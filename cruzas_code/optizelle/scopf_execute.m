function [results, success, raw] = scopf_execute(om, model, mpopt)
%OPF_EXECUTE  Executes the OPF specified by an OPF model object.
%   [RESULTS, SUCCESS, RAW] = OPF_EXECUTE(OM, MPOPT)
%
%   RESULTS are returned with internal indexing, all equipment
%   in-service, etc.
%
%   See also OPF, OPF_SETUP.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
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

%%-----  setup  -----
%% options
dc  = strcmp(upper(mpopt.model), 'DC');
alg = upper(mpopt.opf.ac.solver);
sdp = strcmp(alg, 'SDPOPF');

%% build user-defined costs
om = build_cost_params(om); %for some reason it is necessary


%% get indexing
[vv, ll, nn] = get_idx(om);

if mpopt.verbose > 0
   v = mpver('all');
   fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
end


%%-----  run AC OPF solver  -----
if mpopt.verbose > 0
   fprintf(' -- AC Optimal Power Flow\n');
end

%% if opf.ac.solver not set, choose best available option
if strcmp(alg, 'DEFAULT')
   alg = 'MIPS';             %% MIPS
   mpopt = mpoption(mpopt, 'opf.ac.solver', alg);
end

%% run specific AC OPF solver
switch alg
   case 'MIPS'
      error('MIPS interface not yet implemented with new SCOPF model');
      [results, success, raw] = mipsscopf_solver(om, model, mpopt);
   case 'IPOPT'
      if ~have_fcn('ipopt')
         error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires IPOPT (see http://www.coin-or.org/projects/Ipopt.xml)', alg);
      end
      [results, success, raw] = ipoptscopf_solver(om, model, mpopt);
   case 'OPTIZELLE'
      [results, success, raw] = optizellescopf_solver(om, model, mpopt);
   otherwise
      error('opf_execute: MPOPT.opf.ac.solver = ''%s'' is not a valid AC OPF solver selection', alg);
end

if ~isfield(raw, 'output') || ~isfield(raw.output, 'alg') || isempty(raw.output.alg)
   raw.output.alg = alg;
end