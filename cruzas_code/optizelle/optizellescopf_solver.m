function optizellescopf_solver(mpc, contingencies, mpopt)
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

ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)

% use nominal case to evaluate cost fcn (only pg/qg are relevant)
idx_nom = d.index.getGlobalIndices(mpc, ns, 0);
[VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);

[f, df, d2f] = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), d.om)

end