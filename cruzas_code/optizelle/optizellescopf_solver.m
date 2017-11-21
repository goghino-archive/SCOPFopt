function [results, success, raw] = optizellescopf_solver(om, model, mpopt)
%IPOPTOPF_SOLVER  Solves AC optimal power flow with security constraints using IPOPT.
%
%   [RESULTS, SUCCESS, RAW] = IPOPTSCOPF_SOLVER(OM, MODEL, MPOPT)
%
%   Inputs are an OPF model object, SCOPF model and a MATPOWER options struct.
%
%   Model is a struct with following fields:
%       .cont Containts a list of contingencies
%       .index Contains functions to handle proper indexing of SCOPF variables
%           .getGlobalIndices
%           .getLocalIndicesOPF
%           .getLocalIndicesSCOPF
%           .getREFgens
%           .getPVbuses
%
%   Outputs are a RESULTS struct, SUCCESS flag and RAW output struct.
%
%   The internal x that Optizelle works with has structure
%   [Va1 Vm1 Qg1 Pg_ref1... VaN VmN QgN Pg_refN] [Vm Pg] for all contingency scenarios 1..N
%   with corresponding bounds xmin < x < xmax
%
%   We impose nonlinear equality and inequality constraints g(x) and h(x)
%   with corresponding bounds cl < [g(x); h(x)] < cu
%   and linear constraints l < Ax < u.
%
%   See also OPF, IPOPT.

%   MATPOWER
%   Copyright (c) 2000-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.


%% TODO
% need to work more efficiently with sparse indexing during construction
% of global hessian/jacobian

% how to account for the sparse() leaving out zeros from the sparse
% structure? We want to have exactly same structure across scenarios

%%----- initialization -----
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
   VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
   MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
   QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
   TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
   ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
   deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nn] = get_idx(om);

cont = model.cont;

%% problem dimensions
nb = size(bus, 1);          %% number of buses
ng = size(gen, 1);          %% number of gens
nl = size(branch, 1);       %% number of branches
ns = size(cont, 1);         %% number of scenarios (nominal + ncont)

% get indices of REF gen and of REF/PV buses
[REFgen_idx, nREFgen_idx] = model.index.getREFgens(mpc);
[REFbus_idx,nREFbus_idx] = model.index.getXbuses(mpc,3);%3==REF
[PVbus_idx, nPVbus_idx] = model.index.getXbuses(mpc,2);%2==PV

% indices of local OPF solution vector x = [VA VM PG QG]
[VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);

%% build admittance matrices for nominal case
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% bounds on optimization vars xmin <= x <= xmax
[x0, xmin, xmax] = getv(om); % returns standard OPF form [Va Vm Pg Qg]

% add small pertubation to UB so that we prevent ipopt removing variables
% for which LB=UB, except the Va of the reference bus
tmp = xmax(REFbus_idx);
xmax = xmax + 1e-10;
xmax(REFbus_idx) = tmp;

% replicate bounds for all scenarios and append global limits
xl = xmin([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]); %local variables
xg = xmin([VMopf(PVbus_idx) PGopf(nREFgen_idx)]); %global variables
xmin = [repmat(xl, [ns, 1]); xg];

xl = xmax([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]); %local variables
xg = xmax([VMopf(PVbus_idx) PGopf(nREFgen_idx)]); %global variables
xmax = [repmat(xl, [ns, 1]); xg];

%% try to select an interior initial point based on bounds
if mpopt.opf.init_from_mpc ~= 1
   ll = xmin; uu = xmax;
   ll(xmin == -Inf) = -1e10;               %% replace Inf with numerical proxies
   uu(xmax ==  Inf) =  1e10;
   x0 = (ll + uu) / 2;                     %% set x0 mid-way between bounds
   k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
   x0(k) = xmax(k) - 1;                    %% set just below upper bound
   k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
   x0(k) = xmin(k) + 1;                    %% set just above lower bound
   
   % adjust voltage angles to match reference bus
   Varefs = bus(REFbus_idx, VA) * (pi/180);
   for i = 0:ns-1
      idx = model.index.getGlobalIndices(mpc, ns, i);
      x0(idx(VAscopf)) = Varefs(1);
   end
end

%% find branches with flow limits
il_ = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
il = [1:nl]';               %% we assume every branch has implicit bounds
% TODO insert default limits to branches that
% do not satisfy condition above
nl2 = length(il);           %% number of constrained lines

if size(il_, 1) ~= nl2
   error('Not all branches have specified RATE_A field.');
end

% Define variable bounds.
options.lb = xmin;
options.ub = xmax;

% REVIEW: Could the following (to be modified) by any chance be needed?
% options.cl = [repmat([zeros(2*nb, 1);  -Inf(2*nl2, 1)], [ns, 1]); l];
% options.cu = [repmat([zeros(2*nb, 1); zeros(2*nl2, 1)], [ns, 1]); u+1e10]; %add 1e10 so that ipopt doesn't remove l==u case

% Evaluate cost of nominal case (only Pg/Qg are relevant).
idx_nom = model.index.getGlobalIndices(mpc, ns, 0);

% Pack some additional info to output so that we can verify the solution.
meta.Ybus = Ybus;
meta.Yf = Yf;
meta.Yt = Yt;
meta.lb = options.lb;
meta.ub = options.ub;

%% Optizelle setup.
% Grab the Optizelle library
global Optizelle;
setupOptizelle();

%% Settings for testing
% Toggle on/off if we want to replace infinites with numerical proxies
% (applies only to the variables defined in this file).
replaceInfs = 0;

% Select which optimization we wish to do.
usingUnconstrained = 0;
usingConstrained = 0;
usingEqualityConstrained = 1;
usingInequalityConstrained = 0;

% Select whether we want slack variables or not (only applies to MyEq for
% now).
withSlacks = 0;

%% Further problem settings and computation of (local) minimum.

% Slack variable(s)
s = zeros(ns * 2*nl, 1);

% Bounds on slack variable(s) smin .<= s (no upper bounds)
smin = zeros(size(s));

% Change from infinite bounds to finite bounds xmin .<= x .<= xmax.
if replaceInfs
   xmin(xmin == -Inf) = -1e10;
   xmin(xmin == Inf) = 1e10;
   
   xmax(xmax == -Inf) = -1e10;
   xmax(xmax == Inf) = 1e10;  
end

if withSlacks
   % Append slack variable(s) to initial guess.
   x = [x0; s];
   
   % Append bounds on slack variable(s), considering that there are no 
   % upper bounds on slacks.
   xmin = [xmin; smin];
else
   x = x0;
end

%% Test for Inf, -Inf, and NaN in x.
if find(x == Inf, 1)
   disp('optizellescopf_solver: Inf found in x')
end

if find(x == -Inf, 1)
   disp('optizellescopf_solver: -Inf found in x')
end

if find(isnan(x), 1)
   disp('optizellescopf_solver: NaN found in x')
end

%% Test for Inf, -Inf, and NaN in xmin.
if find(xmin == Inf, 1)
   disp('optizellescopf_solver: Inf found in xmin')
end

if find(xmin == -Inf, 1)
   disp('optizellescopf_solver: -Inf found in xmin')
end

if find(isnan(xmin), 1)
   disp('optizellescopf_solver: NaN found in xmin')
end

%% Test for Inf, -Inf, and NaN in xmax.
if find(xmax == Inf, 1)
   disp('optizellescopf_solver: Inf found in xmax')
end

if find(xmax == -Inf, 1)
   disp('optizellescopf_solver: -Inf found in xmax')
end

if find(isnan(xmax), 1)
   disp('optizellescopf_solver: NaN found in xmax')
end

% Number of equality constraints in g_new(x) = [g(x), h(x) - z].
% g(x) = 0, h(x) >= 0, z >= 0.
%TODO: inequality constraints are formulated as h(x) <= 0, see opf_consfcn()

if withSlacks
   NEQ = 2*nb + 2*nl;
else
   NEQ = 2*nb;
end

% Number of inequality constraints (xmin .<= x_with_s, x_no_s .<= xmax),
% where xmin, x_with_s, x_no_s, and xmax are (column) vectors.
if withSlacks
   NINEQ = 2*length(x0) + length(s);
else
   NINEQ = 2*nl;
end

myauxdata.idx_nom = idx_nom;
myauxdata.model = model;
myauxdata.om = om;
myauxdata.mpc = mpc;
myauxdata.il = il;
myauxdata.mpopt = mpopt;
myauxdata.xmin = xmin;
myauxdata.xmax = xmax;
myauxdata.NEQ = NEQ;
myauxdata.NINEQ = NINEQ;
myauxdata.lenx_no_s = length(x0); % length of x without slack variables.
myauxdata.withSlacks = withSlacks;

% Allocate memory for the equality multiplier
y = zeros(ns*NEQ, 1);

% Allocate memory for the inequality multiplier
z = zeros(NINEQ, 1);

if usingUnconstrained
   disp('Using Unconstrained...')
elseif usingConstrained
   disp('Using Constrained...')
elseif usingEqualityConstrained
   disp('Using EqualityConstrained...')
elseif usingInequalityConstrained 
   disp('Using InequalityConstrained...')
end

if withSlacks
   disp('Using slack variables...')
end

% No need to check if more than one selected since using if...elseif statements
% below.

if usingUnconstrained

   % Create an optimization state
   state = Optizelle.Unconstrained.State.t(Optizelle.Rm, x);
   
   % Create a bundle of functions
   fns = Optizelle.Unconstrained.Functions.t;
   
elseif usingConstrained
   
   % Create an optimization state
   state = Optizelle.Constrained.State.t( ...
      Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,x,y,z);
   
   % Create a bundle of functions
   fns = Optizelle.Constrained.Functions.t;
   
elseif usingEqualityConstrained
   
   % Create an optimization state
   state = Optizelle.EqualityConstrained.State.t( ...
      Optizelle.Rm,Optizelle.Rm,x,y);
   
   % Create a bundle of functions
   fns = Optizelle.EqualityConstrained.Functions.t;
   
elseif usingInequalityConstrained
   
   state=Optizelle.InequalityConstrained.State.t( ...
      Optizelle.Rm,Optizelle.Rm,x,z);
   
   % Create a bundle of functions
   fns=Optizelle.InequalityConstrained.Functions.t;
end

% Define bundle of functions.
fns.f = MyObj(myauxdata);

if usingConstrained || usingEqualityConstrained
   fns.g = MyEq(myauxdata);
end

if usingConstrained || usingInequalityConstrained
   fns.h = MyIneq(myauxdata);
end

% Solve the optimization problem
% tic
if usingUnconstrained

    state = Optizelle.Unconstrained.Algorithms.getMin( ...
        Optizelle.Rm, Optizelle.Messaging.stdout, ...
        fns,state);
   
elseif usingConstrained
   
   state = Optizelle.Constrained.Algorithms.getMin( ...
      Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout, ...
      fns,state);
   
elseif usingEqualityConstrained
   
   state = Optizelle.EqualityConstrained.Algorithms.getMin( ...
      Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout, ...
      fns,state);
   
elseif usingInequalityConstrained
   
   state = Optizelle.InequalityConstrained.Algorithms.getMin( ...
      Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);
end

% toc

% Print out the reason for convergence
fprintf('The algorithm converged due to: %s\n', ...
   Optizelle.OptimizationStop.to_string(state.opt_stop));

% TODO: Figure out exactly what to return in results, success, and raw.
% TODO: in raw, make sure you put the time taken by the algorithm.
% Check line 389 of ipoptscopf_solver.
raw = struct('meta', meta);
results = struct('x', state.x);
success = 1;

% raw = struct('info', info.status, 'meta', meta, 'numIter', info.iter, 'overallAlgorithm', info.cpu);
% results = struct('f', f, 'x', x);
end
