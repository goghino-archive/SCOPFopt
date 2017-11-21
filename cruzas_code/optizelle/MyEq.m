%% Define equality constraints.
function self = MyEq(myauxdata)

% y=g(x)
self.eval = @(x) constraints(x, myauxdata);

% y=g'(x)dx
self.p = @(x,dx) jacobvec(x, dx, myauxdata);

% xhat=g'(x)*dy
self.ps = @(x,dy) jacobvec(x, dy, myauxdata);

% xhat=(g''(x)dx)*dy
self.pps = @(x,dx,dy) hessian(x, myauxdata, dy) * dx;
% self.pps = @(x,dx,dy) 0;

% Helper functions.
   function constr = constraints(x, myauxdata)
      % Extract data.
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      model = myauxdata.model;
      mpopt = myauxdata.mpopt;
      il = myauxdata.il;
      lenx_no_s = myauxdata.lenx_no_s; % length of x without slack variables.
      NEQ = myauxdata.NEQ; % number of equality constraints.
      withSlacks = myauxdata.withSlacks;
      
      nb = size(mpc.bus, 1);     %% number of buses
      nl = size(mpc.branch, 1);  %% number of branches
      ns = size(model.cont, 1);  %% number of scenarios (nominal + ncont)
      
      constr = zeros(ns*(NEQ), 1);
      
      [VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
      
      for i = 0:ns-1
         cont = model.cont(i+1);
         idx = model.index.getGlobalIndices(mpc, ns, i);
         
         [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
         [hn_local, gn_local] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), om, Ybus, Yf, Yt, mpopt, il);
         
         if withSlacks
            % Extract slack variable(s) s from x.
            % s = [s0, s1, ... , sNS]
            s = x(lenx_no_s + i*2*nl + (1:2*nl));
         
            % Since h(x) <= 0 in Matpower, -h(x) - s = 0, s >= 0
            % as required by Optizelle.
            constr(i*(NEQ) + (1:NEQ)) = [gn_local; -hn_local - s];
         else
            constr(i*(NEQ) + (1:NEQ)) = gn_local;
         end
      end
      
      
      %% Test for Inf, -Inf, and NaN in constr.
      if find(constr == Inf, 1)
         disp('MyEq: Inf found in constr')
      end
      
      if find(constr == -Inf, 1)
         disp('MyEq: -Inf found in constr')
      end
      
      if find(isnan(constr), 1)
         disp('MyEq: NaN found in constr')
      end
   end

   function jvec = jacobvec(x, d, myauxdata)
      % Extract data.
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      model = myauxdata.model;
      mpopt = myauxdata.mpopt;
      il = myauxdata.il;
      lenx_no_s = myauxdata.lenx_no_s; % length of x without slack variables.
      NEQ = myauxdata.NEQ; % number of equality constraints.
      withSlacks = myauxdata.withSlacks;
      
      
      nb = size(mpc.bus, 1);     %% number of buses
      nl = size(mpc.branch, 1);  %% number of branches
      ns = size(model.cont, 1);     %% number of scenarios (nominal + ncont)
      
      % get indices of REF gen and PV bus
      [REFgen_idx, nREFgen_idx] = model.index.getREFgens(mpc);
      [PVbus_idx, nPVbus_idx] = model.index.getXbuses(mpc,2);%2==PV
      
      [VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
      [VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);
      
      % -Id. Placed in the lower-right "corner" of each scenario's
      % Jacobian matrix.
      neg_identity = -sparse(eye(2*nl, 2*nl));
      
      J = sparse(ns*NEQ, length(x));
      
      for i = 0:ns-1
         idx = model.index.getGlobalIndices(mpc, ns, i);
         
         cont = model.cont(i+1);
         [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
         [hn, gn, dhn, dgn] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), om, Ybus, Yf, Yt, mpopt, il);
         
         % Transpose since opf_consfcn transposed solutions.
         % Take negative since Matpower requires h(x) <= 0 but
         % Optizelle requires h(x) >= 0.
         dhn = -dhn';
         dgn = dgn';
         
         if withSlacks
            %jacobian wrt local variables
            J(i*NEQ + (1:NEQ), idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)])) = [dgn(:,[VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]);...
               dhn(:,[VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)])];
            %jacobian wrt global variables
            J(i*NEQ + (1:NEQ), idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)])) = [dgn(:, [VMopf(PVbus_idx) PGopf(nREFgen_idx)]);...
               dhn(:, [VMopf(PVbus_idx) PGopf(nREFgen_idx)])];
            
         else
            %jacobian wrt local variables
            J(i*NEQ + (1:NEQ), idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)])) = dgn(:,[VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]);
            
            %jacobian wrt global variables
            J(i*NEQ + (1:NEQ), idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)])) = dgn(:, [VMopf(PVbus_idx) PGopf(nREFgen_idx)]);
         end
         
         % Set corner of Jacobian of this scenario to -Id.
         % Number of rows in dgn is 2*nb; in dhn it is 2*nl.
         % Number of columns in dgn, and dhn, is lenx_no_s.
         % Structure of Jacobian in scenario i is:
         % J = [dgn, 0; dhn, -Id]. Hence why we use following row and
         % column indices.
         % Note: -Id has dimensions 2*nl by 2*nl.
         %
         % Considering the system of slack variables, we have
         % [g(x); h(x) - s] for the equality constraints.
         % As we can see, we have 0s in the upper right corner of J since the
         % partial derivative of the constraints g(x) w.r.t. the slack variables
         % is 0, since g(x) does not at all depend on s.
         % Similarly, we have -Id in the specified (lower-right) corner
         % given the constraints w.r.t. the slack variables are
         % h(x) - s: taking the partial derivative of h(x) - s w.r.t. the slack
         % variables yields the -Id matrix.
         if withSlacks
            J(i*NEQ + 2*nb + (1:2*nl), lenx_no_s + i*2*nl + (1:2*nl)) = neg_identity;
         end
      end
      
      dType = 0;
      if (size(d, 1) == size(x, 1))  % case: d == dx
         dType = 'dx';
         jvec = J * d;
      elseif (size(d, 1) == ns*NEQ) % case: d == dy
         jvec = J' * d;
         dType = 'dy';
      end
      
      %% Test for Inf, -Inf, and NaN in d.
      if find(d == Inf)
         fprintf('MyEq: Inf found in %s\n', dType);
      end
      
      if find(d == -Inf)
         fprintf('MyEq: -Inf found in %s\n', dType);
      end
      
      if find(isnan(d))
         fprintf('MyEq: NaN found in %s\n', dType);
      end
      
      %% Test for Inf, -Inf, and NaN in jvec.
      if find(jvec == Inf, 1)
         disp('MyEq: Inf found in jvec')
      end
      
      if find(jvec == -Inf, 1)
         disp('MyEq: -Inf found in jvec')
      end
      
      if find(isnan(jvec), 1)
         disp('MyEq: NaN found in jvec')
      end
   end

   function H = hessian(x, myauxdata, dy)
      % Extract data.
      il = myauxdata.il;
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      model = myauxdata.model;
      mpopt = myauxdata.mpopt;
      NEQ = myauxdata.NEQ; % number of equality constraints.
      NINEQ = myauxdata.NINEQ;
      withSlacks = myauxdata.withSlacks;
      
      nb = size(mpc.bus, 1);          %% number of buses
      nl = size(mpc.branch, 1);       %% number of branches
      ns = size(model.cont, 1);           %% number of scenarios (nominal + ncont)
      
      H = sparse(size(x,1), size(x,1));
      
      % get indices of REF gen and PV bus
      [REFgen_idx, nREFgen_idx] = model.index.getREFgens(mpc);
      [PVbus_idx, nPVbus_idx] = model.index.getXbuses(mpc,2);%2==PV
      
      [VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
      [VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);
      
      % REVIEW: is this correct?
      sigma = 0;
      
      for i = 0:ns-1
         %compute local indices and its parts
         idx = model.index.getGlobalIndices(mpc, ns, i);
         
         cont = model.cont(i+1);
         [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
         
         % c(x) = [gn0; hn0 - s0; gn1; hn1 - s1; ... ; gnNs; hnNs - sNs]
         % (Ns is the number of scenarios)
         % where gn corresponds to equality constraints, hn corresponds to
         % inequality constraints.
         %
         % The lagrange multipliers in the dy will be composed accordingly:
         % dy = [lamEq0; lamIneq0; lamEq1; lamIneq1; ... ; lamEqNs; lamIneqNs]
         % Hence why we need to extract the lagrange multipliers as
         % follows.
         lam.eqnonlin = dy(i*NEQ + (1:2*nb), 1);
         
         if withSlacks
            lam.ineqnonlin = dy(i*NEQ + 2*nb + (1:2*nl), 1);
         else
            lam.ineqnonlin = zeros(NINEQ, 1);
         end
         
         % Take negative since Matpower requires h(x) <= 0 but
         % Optizelle requires h(x) >= 0.
         lam.ineqnonlin = -lam.ineqnonlin;
            
         H_local = opf_hessfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), lam, sigma, om, Ybus, Yf, Yt, mpopt, il);
         
         % H_ll (PG_ref relevant only in nominal case, added to global part)
         H(idx([VAscopf VMscopf(nPVbus_idx) QGscopf]), idx([VAscopf VMscopf(nPVbus_idx) QGscopf])) =...
            H_local([VAopf VMopf(nPVbus_idx) QGopf], [VAopf VMopf(nPVbus_idx) QGopf]);
         
         % H_lg and H_gl (PG parts are implicitly zero, could leave them out)
         H(idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)]), idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)])) = ...
            H_local([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)], [VMopf(PVbus_idx) PGopf(nREFgen_idx)]);
         H(idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)]), idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)])) = ...
            H_local([VMopf(PVbus_idx) PGopf(nREFgen_idx)], [VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]);
         
         % H_gg hessian w.r.t global variables (and PG_ref_0)
         if i == 0
            % H_pg at non-reference gens, these are global variables
            H(idx([PGscopf(nREFgen_idx)]), idx([PGscopf(nREFgen_idx)])) = ...
               H_local([PGopf(nREFgen_idx)], [PGopf(nREFgen_idx)]);
            
            % H_pgref is local variable for nominal scenario, but used in f()
            H(idx([PGscopf(REFgen_idx)]), idx([PGscopf(REFgen_idx)])) = ...
               H_local([PGopf(REFgen_idx)], [PGopf(REFgen_idx)]);
         end
         
         %each scenario contributes to hessian w.r.t global VM variables at PV buses
         H(idx([VMscopf(PVbus_idx)]), idx([VMscopf(PVbus_idx)])) = ...
            H(idx([VMscopf(PVbus_idx)]), idx([VMscopf(PVbus_idx)])) + ...
            H_local([VMopf(PVbus_idx)], [VMopf(PVbus_idx)]);
      end
      
      %% Test for Inf, -Inf, and NaN in H.
      if find(H == Inf, 1)
         disp('MyEq: Inf found in H')
      end
      
      if find(H == -Inf, 1)
         disp('MyEq: -Inf found in H')
      end
      
      if find(isnan(H), 1)
         disp('MyEq: NaN found in H')
      end
   end
end