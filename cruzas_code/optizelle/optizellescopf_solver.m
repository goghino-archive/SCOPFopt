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

%% build local connectivity matrices
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
Cl = Cf + Ct;                                   %% for each line - from & to
Cb_nominal = Cl' * Cl + speye(nb);              %% for each bus - contains adjacent buses
Cl2_nominal = Cl(il, :);                        %% branches with active flow limit
Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng); %%locations where each gen. resides

% Define variable bounds.
options.lb = xmin;
options.ub = xmax;

idx_nom = model.index.getGlobalIndices(mpc, ns, 0); %evaluate cost of nominal case (only Pg/Qg are relevant)

%pack some additional info to output so that we can verify the solution
meta.Ybus = Ybus;
meta.Yf = Yf;
meta.Yt = Yt;

% Grab the Optizelle library
global Optizelle;
setupOptizelle();

% Initial guess.
x = x0;
% TODO: must have equality and inequality multipliers when we consider
% constraints.

myauxdata.idx_nom = idx_nom;
myauxdata.model = model;
myauxdata.om = om;
myauxdata.mpc = mpc;

% Create an optimization state
state = Optizelle.Unconstrained.State.t(Optizelle.Rm, x);

% Create a bundle of functions
fns = Optizelle.Unconstrained.Functions.t;
fns.f = MyObj(myauxdata);

% Solve the optimization problem
state = Optizelle.Unconstrained.Algorithms.getMin( ...
   Optizelle.Rm, Optizelle.Messaging.stdout, ...
   fns, state);

% Print out the reason for convergence
fprintf('The algorithm converged due to: %s\n', ...
   Optizelle.OptimizationStop.to_string(state.opt_stop));

% Pack some additional info to output so that we can verify the solution.
meta.Ybus = Ybus;
meta.Yf = Yf;
meta.Yt = Yt;
meta.lb = options.lb;
meta.ub = options.ub;

% TODO: Figure out exactly what to return in results, success, and raw.
% TODO: in raw, make sure you put the time taken by the algorithm.
% Check line 389 of ipoptscopf_solver.
raw = struct('meta', meta);
results = struct('x', state.x);
success = 1;
end

%% Define objective function.
function self = MyObj(myauxdata)
% Evaluation
self.eval = @(x) myFEval(x, myauxdata);

% Gradient
self.grad = @(x) myDfEval(x, myauxdata);

% Hessian-vector product
self.hessvec = @(x, dx) myD2fEval(x, myauxdata) * dx;


% Helper functions.
   function f = myFEval(x, myauxdata)
      % Extract data.
      idx_nom = myauxdata.idx_nom;
      model = myauxdata.model;
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      
      % Indices of local OPF solution vector x = [VA VM PG QG].
      [VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
      
      f = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), om);
   end

   function grad = myDfEval(x, myauxdata)
      % Extract data.
      idx_nom = myauxdata.idx_nom;
      model = myauxdata.model;
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      
      % Indices of local OPF solution vector x = [VA VM PG QG].
      [VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
      [VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);
      
      [f, df] = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), om);
      
      % Nonzero only nominal case Pg.
      grad = zeros(size(x,1),1);
      grad(idx_nom(PGscopf)) = df(PGopf);
   end

   function hess = myD2fEval(x, myauxdata)
      % Extract data.
      idx_nom = myauxdata.idx_nom;
      model = myauxdata.model;
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      
      % Indices of local OPF solution vector x = [VA VM PG QG].
      [VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
      [VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);
      
      [f, df, d2f] = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), om);
      
      % Nonzero only nominal case Pg.
      hess = sparse(size(x,1),size(x,1));
      hess(idx_nom(PGscopf), idx_nom(PGscopf)) = d2f(PGopf, PGopf);
   end


end

%% Define equality constraints.
function self = MyEq()
%
%     % y=g(x)
%     self.eval = @(x) sum(x.^2) - 40;
%
%     % y=g'(x)dx
%     self.p = @(x,dx) [2*x(1), 2*x(2), 2*x(3), 2*x(4)] * dx;
%
%     % xhat=g'(x)*dy
%     self.ps = @(x,dy)  [2*x(1); 2*x(2); 2*x(3); 2*x(4)] .* dy;
%
%     % xhat=(g''(x)dx)*dy
%     self.pps = @(x,dx,dy) [2*dx(1); 2*dx(2); 2*dx(3); 2*dx(4)] .* dy;
end


function constr = constraints(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
NCONSTR = 2*nb + 2*nl;

constr = zeros(ns*(NCONSTR), 1);

[VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

for i = 0:ns-1
   cont = d.cont(i+1);
   idx = d.index.getGlobalIndices(mpc, ns, i);
   [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
   [hn_local, gn_local] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), d.om, Ybus, Yf, Yt, d.mpopt, d.il);
   constr(i*(NCONSTR) + (1:NCONSTR)) = [gn_local; hn_local];
end

if ~isempty(d.A)
   constr = [constr; d.A*x]; %append linear constraints
end
end

function J = jacobian(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
NCONSTR = 2*nb + 2*nl;          %% number of constraints (eq + ineq)

J = sparse(ns*(NCONSTR), size(x,1));

% get indices of REF gen and PV bus
[REFgen_idx, nREFgen_idx] = d.index.getREFgens(mpc);
[PVbus_idx, nPVbus_idx] = d.index.getXbuses(mpc,2);%2==PV

[VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

for i = 0:ns-1
   %compute local indices
   idx = d.index.getGlobalIndices(mpc, ns, i);
   
   cont = d.cont(i+1);
   [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
   [hn, gn, dhn, dgn] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), d.om, Ybus, Yf, Yt, d.mpopt, d.il);
   dgn = dgn';
   dhn = dhn';
   
   %jacobian wrt local variables
   J(i*NCONSTR + (1:NCONSTR), idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)])) = [dgn(:,[VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]);...
      dhn(:,[VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)])];
   %jacobian wrt global variables
   J(i*NCONSTR + (1:NCONSTR), idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)])) = [dgn(:, [VMopf(PVbus_idx) PGopf(nREFgen_idx)]);...
      dhn(:, [VMopf(PVbus_idx) PGopf(nREFgen_idx)])];
end
J = [J; d.A]; %append Jacobian of linear constraints
end