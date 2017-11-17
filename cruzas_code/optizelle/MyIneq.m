%% Define inequalities.
function self = MyIneq(myauxdata)
% z=h(x)
self.eval = @(x) constraints(x, myauxdata);

% z=h'(x)dx
self.p = @(x,dx) jacobvec(x, dx);

% xhat=h'(x)*dz
% self.ps = @(x,dz) jacobvec(x, dz);
self.ps = @(x,dz) 0;

% xhat=(h''(x)dx)*dz
self.pps = @(x,dx,dz) sparse(length(x),length(x)); 


% Helper functions.
   function constr = constraints(x, myauxdata)
      % Append constraints for xmin <= x <= xmax.
      % -> x - xmin >= 0
      % -> xmax - x >= 0
      constr = [x - myauxdata.xmin;
                myauxdata.xmax - x];
   end

   function jvec = jacobvec(x, d)
      x_size = size(x)
      d_size = size(d)
      
     % Append Jacobian for x - xmin >=0.
     J = sparse(eye(length(x)));

     % Append Jacobian for xmax - x >= 0.
     J = [J; -sparse(eye(length(x)))];
     
     J_size = size(J)
     
     
      if (size(d, 1) == size(x, 1))  % case: d == dx
         jvec = J * d;
      elseif (size(d, 1) == myauxdata.NINEQ) % case: d == dz
         jvec = J' * d;
      end
     
     jvec_size = size(jvec)
   end
end