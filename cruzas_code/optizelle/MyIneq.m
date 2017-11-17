%% Define inequalities.
function self = MyIneq(myauxdata)
% z=h(x)
self.eval = @(x) constraints(x, myauxdata);

% z=h'(x)dx
self.p = @(x,dx) jacobvec(x, dx);
% self.p = @(x,dz) 0;

% xhat=h'(x)*dz
self.ps = @(x,dz) jacobvec(x, dz);
% self.ps = @(x,dz) 0;

% xhat=(h''(x)dx)*dz
self.pps = @(x,dx,dz) sparse(length(x),length(x));
% self.pps = @(x,dx,dz) 0;


% Helper functions.
   function constr = constraints(x, myauxdata)
      % Append constraints for xmin <= x <= xmax.
      xmin = myauxdata.xmin;
      xmax = myauxdata.xmax;
      
%       xmin(xmin == -Inf) = 1e-12;
%       xmin(xmin == Inf) = 1e20;
%       
%       xmax(xmax == -Inf) = 1e-12;
%       xmax(xmax == Inf) = 1e20;
      
      % -> x - xmin >= 0
      % -> xmax - x >= 0
      constr = [x - xmin;
                xmax - x];
      %       constr = [x - (1e-6 + zeros(length(x),1));
      %                 1e4*ones(length(x),1) - x];
   end

   function jvec = jacobvec(x, d)
      % Append Jacobian for x - xmin >=0.
      J = sparse(eye(length(x)));
      
      % Append Jacobian for xmax - x >= 0.
      J = [J; -sparse(eye(length(x)))];
      
      if (size(d, 1) == size(x, 1))  % case: d == dx
         %          disp('d == dx')
         jvec = J * d;
      elseif (size(d, 1) == myauxdata.NINEQ) % case: d == dz
         %          disp('d == dz')
         jvec = J' * d;
      end
      
      %      jvec_size = size(jvec)
   end
end