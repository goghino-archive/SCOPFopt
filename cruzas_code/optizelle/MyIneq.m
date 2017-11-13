%% Define inequalities.
function self = MyIneq(myauxdata)
% z=h(x)
self.eval = @(x) constraints(x, myauxdata);

% z=h'(x)dx
% self.p = @(x,dx) jacobian(x) * dx;
self.p = @(x,dx) bs(jacobian(x), dx);

% xhat=h'(x)*dz
% self.ps = @(x,dz) jacobian(x)' * dz;
self.ps = @(x,dz) bs(jacobian(x)', dz);

% xhat=(h''(x)dx)*dz
self.pps = @(x,dx,dz) [ 0. ]; 

   function res = bs(A, x)
      A_size = size(A)
      x_size = size(x)
      res = 0;
   end

% Helper functions.
   function constr = constraints(x, myauxdata)
      % Append constraints for xmin <= x <= xmax.
      constr = [ x >= myauxdata.xmin;
                myauxdata.xmax >= x];
   end

   function J = jacobian(x)
     % Append Jacobian for x >= xmin.
     J = sparse(eye(length(x)));

     % Append Jacobian for x <= xmax.
     J = [J; -sparse(eye(length(x)))];
   end
end