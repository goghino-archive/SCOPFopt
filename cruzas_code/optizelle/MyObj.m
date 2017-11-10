%% Define objective function.
function self = MyObj(myauxdata)
% Evaluation
self.eval = @(x) objective(x, myauxdata);

% Gradient
self.grad = @(x) gradient(x, myauxdata);

% Hessian-vector product
self.hessvec = @(x, dx) hessian(x, myauxdata) * dx;

% Helper functions.
   function f = objective(x, myauxdata)
      % Extract data.
      idx_nom = myauxdata.idx_nom;
      model = myauxdata.model;
      om = myauxdata.om;
      mpc = myauxdata.mpc;
      
      % Indices of local OPF solution vector x = [VA VM PG QG].
      [VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
      
      f = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), om);
   end

   function grad = gradient(x, myauxdata)
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

   function H = hessian(x, myauxdata)
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
      H = sparse(size(x,1),size(x,1));
      H(idx_nom(PGscopf), idx_nom(PGscopf)) = d2f(PGopf, PGopf);
   end
end