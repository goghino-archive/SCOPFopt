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
      NINEQ = myauxdata.NINEQ;
      
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
         s = x(lenx_no_s+1 : end);
         
         constr(i*(NEQ + NINEQ) + (1:NEQ)) = gn_local;
         constr(i*(NEQ + NINEQ) + (NEQ+1:NEQ+NINEQ)) = hn_local - s;
      end
   end

   function J = jacobian(x, myauxdata)
      % Extract data.
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      model = myauxdata.model;
      mpopt = myauxdata.mpopt;
      il = myauxdata.il;
      NEQ = myauxdata.NEQ; % number of equality constraints.
      
      nb = size(mpc.bus, 1);     %% number of buses
      ns = size(model.cont, 1);     %% number of scenarios (nominal + ncont)
      
      J = sparse(ns*(NEQ), size(x,1));
      
      % get indices of REF gen and PV bus
      [REFgen_idx, nREFgen_idx] = model.index.getREFgens(mpc);
      [PVbus_idx, nPVbus_idx] = model.index.getXbuses(mpc,2);%2==PV
      
      [VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
      [VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);
      
      for i = 0:ns-1
         idx = model.index.getGlobalIndices(mpc, ns, i);
         
         cont = model.cont(i+1);
         [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
         [hn, gn, dhn, dgn] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), om, Ybus, Yf, Yt, mpopt, il);
         
         % Transpose since opf_consfcn returns dgn'.
         dgn = dgn';
         
         %jacobian wrt local variables
         J(i*NEQ + (1:NEQ), ...
            idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)])) =...
            dgn(:,[VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]);
         
         %jacobian wrt global variables
         J(i*NEQ + (1:NEQ), ...
            idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)])) =...
            dgn(:, [VMopf(PVbus_idx) PGopf(nREFgen_idx)]);
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
      
      lam.ineqnonlin = zeros(2*nl, 1);
      sigma = 0;
      
      for i = 0:ns-1
         %compute local indices and its parts
         idx = model.index.getGlobalIndices(mpc, ns, i);
         
         cont = model.cont(i+1);
         [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
         
         lam.eqnonlin = dy(i*NEQ + (1:NEQ), 1);
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