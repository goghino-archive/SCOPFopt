%% Define equality constraints.
function self = MyEq(myauxdata)

% y=g(x)
self.eval = @(x) constraints(x, myauxdata);

% y=g'(x)dx
self.p = @(x,dx) jacobian(x, myauxdata) * dx;

% xhat=g'(x)*dy
self.ps = @(x,dy) jacobian(x, myauxdata)' * dy;

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
         
         % Extract slack variable(s) s from x.
         s = x(lenx_no_s + (1:2*nl));
         
         constr(i*(NEQ) + (1:NEQ)) = [gn_local; hn_local - s];
      end
   end

   function J = jacobian(x, myauxdata)
      % Extract data.
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      model = myauxdata.model;
      mpopt = myauxdata.mpopt;
      il = myauxdata.il;
      lenx_no_s = myauxdata.lenx_no_s; % length of x without slack variables.
      NEQ = myauxdata.NEQ; % number of equality constraints.
      
      
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
      
      for i = 0:ns-1
         idx = model.index.getGlobalIndices(mpc, ns, i);
         
         cont = model.cont(i+1);
         [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
         [hn, gn, dhn, dgn] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), om, Ybus, Yf, Yt, mpopt, il);
         
         % Transpose since opf_consfcn transposed solutions.
         dhn = dhn';
         dgn = dgn';
         
         %jacobian wrt local variables
         J(i*NEQ + (1:NEQ), idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)])) = [dgn(:,[VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]);...
            dhn(:,[VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)])];
         %jacobian wrt global variables
         J(i*NEQ + (1:NEQ), idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)])) = [dgn(:, [VMopf(PVbus_idx) PGopf(nREFgen_idx)]);...
            dhn(:, [VMopf(PVbus_idx) PGopf(nREFgen_idx)])];
         
         % Set corner of Jacobian of this scenario to -Id.
         % Number of rows in dgn is 2*nb; in dhn it is 2*nl.
         % Number of columns in dgn, and dhn, is lenx_no_s. 
         % Structure of Jacobian in scenario i is:
         % J = [dgn, 0s ; dhn -Id]. Hence why we use following row and
         % column indices.
         % Note: -Id has dimensions 2*nl by 2*nl.
         J(i*NEQ + 2*nb + (1:2*nl), i*NEQ + lenx_no_s + (1:2*nl)) = neg_identity;
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
         lam.ineqnonlin = dy(i*NEQ + 2*nb + (1:2*nl), 1);
         
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
   end
end