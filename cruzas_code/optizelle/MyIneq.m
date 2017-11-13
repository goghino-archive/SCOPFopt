%% Define inequalities.
function self = MyIneq(myauxdata)
% z=h(x)
self.eval = @(x) constraints(x, myauxdata);

% z=h'(x)dx
self.p = @(x,dx) jacobian(x, myauxdata) * dx;

% xhat=h'(x)*dz
self.ps = @(x,dz) jacobian(x, myauxdata)' * dz;

% xhat=(h''(x)dx)*dz
self.pps = @(x,dx,dz) 0;

   % Function for checking dimensions.
   function res = bs(x, A, dx)
      x_size = size(x)
      A_size = size(A)
      dx_size = size(dx)
      Adx_size = size(A*dx)
      res = 0;
   end

% Helper functions.
   function constr = constraints(x, myauxdata)
      % Append constraints for xmin <= x <= xmax.
      constr = [ x >= myauxdata.xmin;
                myauxdata.xmax >= x];
   end

   function J = jacobian(x, myauxdata)
     % Append Jacobian for x >= xmin.
     J = sparse(eye(length(x)));

     % Append Jacobian for x <= xmax.
     J = [J; -sparse(eye(length(x)))];
   end
end