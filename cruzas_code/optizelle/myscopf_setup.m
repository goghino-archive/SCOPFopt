function om = myscopf_setup(mpc)
%OPF  Constructs an OPF model object from a MATPOWER case struct.
%   OM = OPF_SETUP(MPC, MPOPT)
%
%   Assumes that MPC is a MATPOWER case struct with internal indexing,
%   all equipment in-service, etc.
%
%   See also OPF, EXT2INT, OPF_EXECUTE.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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

%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
if isfield(mpc, 'A')
  error('User constraints not supported in SCOPF');
end

if isfield(mpc, 'N')
    error('general cost vars not supported in SCOPF');
end

%% create (read-only) copies of individual fields for convenience
[baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);

%% warn if there is more than one reference bus
refs = find(bus(:, BUS_TYPE) == REF);
if length(refs) > 1 && mpopt.verbose > 0
  errstr = ['\nopf_setup: Warning: Multiple reference buses.\n', ...
              '           For a system with islands, a reference bus in each island\n', ...
              '           may help convergence, but in a fully connected system such\n', ...
              '           a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end

%% set up initial variables and bounds
Va   = bus(:, VA) * (pi/180);
Vm   = bus(:, VM);
Pg   = gen(:, PG) / baseMVA;
Qg   = gen(:, QG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
Qmin = gen(:, QMIN) / baseMVA;
Qmax = gen(:, QMAX) / baseMVA;


%% voltage angle reference constraints
Vau = Inf(nb, 1);
Val = -Vau;
Vau(refs) = Va(refs);
Val(refs) = Va(refs);

if any(gencost(:, MODEL) ~= POLYNOMIAL & gencost(:, MODEL) ~= PW_LINEAR)
    error('opf_setup: some generator cost rows have invalid MODEL value');
end

%% construct OPF model object
om = opf_model(mpc);

om = add_vars(om, 'Va', nb, Va, Val, Vau); %adds OPF variables for voltage angles and set initial values
om = add_vars(om, 'Vm', nb, Vm, bus(:, VMIN), bus(:, VMAX)); %voltage magnitude
om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax); %generator real power
om = add_vars(om, 'Qg', ng, Qg, Qmin, Qmax); %generator imaginary power
om = add_constraints(om, 'Pmis', nb, 'nonlinear');
om = add_constraints(om, 'Qmis', nb, 'nonlinear');%real and reactive power balance (eq. constr.)
om = add_constraints(om, 'Sf', nl, 'nonlinear');
om = add_constraints(om, 'St', nl, 'nonlinear');%branch flow limits (ineq. constr.)
