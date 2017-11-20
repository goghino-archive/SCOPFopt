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
      
      % -> x - xmin >= 0
      % -> xmax - x >= 0
      constr = [x - xmin;
                xmax - x];
%             constr = [x - (1e-6 + zeros(length(x),1));
%                       1e4*ones(length(x),1) - x];

      %% Test for Inf, -Inf, and NaN in constr.
      if find(constr == Inf, 1)
         disp('MyIneq: Inf found in constr')
      end
      
      if find(constr == -Inf, 1)
         disp('MyIneq: -Inf found in constr')
      end
      
      if find(isnan(constr), 1)
         disp('MyIneq: NaN found in constr')
      end
   end

   function jvec = jacobvec(x, d)
      % Append Jacobian for x - xmin >=0.
      J = sparse(eye(length(x)));
      
      % Append Jacobian for xmax - x >= 0.
      J = [J; -sparse(eye(length(x)))];
      
      % Replace infinite bounds in d with finite bounds.
%       d(d == -Inf) = -1e10;
%       d(d == Inf) = 1e10;
% 
%       d(d == -Inf) = -1e10;
%       d(d == Inf) = 1e10;
      
      dType = 0;
      if (size(d, 1) == size(x, 1))  % case: d == dx
         jvec = J * d;
         dType = 'dx';
      elseif (size(d, 1) == myauxdata.NINEQ) % case: d == dz
         dType = 'dz';
         jvec = J' * d;
      end
      
%       jvec(isnan(jvec)) = 0;
      
      %% Test for Inf, -Inf, and NaN in d.
      if find(d == Inf)
         fprintf('MyIneq: Inf found in %s\n', dType);
      end
      
      if find(d == -Inf)
         fprintf('MyIneq: -Inf found in %s\n', dType);
      end
      
      if find(isnan(d))
         fprintf('MyIneq: NaN found in %s\n', dType);
      end
      
      %% Test for Inf, -Inf, and NaN in jvec.
      if find(jvec == Inf, 1)
         disp('MyIneq: Inf found in jvec')
      end
      
      if find(jvec == -Inf, 1)
         disp('MyIneq: -Inf found in jvec')
      end
      
      if find(isnan(jvec), 1)
         disp('MyIneq: NaN found in jvec')
      end
   end
end