%% Define inequalities.
function self = MyIneq(myauxdata)
% z=h(x)
self.eval = @(x) constraints(x, myauxdata);

% z=h'(x)dx
self.p = @(x,dx) jacobvec(x, dx, myauxdata);
% self.p = @(x,dz) 0;

% xhat=h'(x)*dz
self.ps = @(x,dz) jacobvec(x, dz, myauxdata);
% self.ps = @(x,dz) 0;

% xhat=(h''(x)dx)*dz
self.pps = @(x,dx,dz) sparse(length(x),length(x));
% self.pps = @(x,dx,dz) 0;


% Helper functions.
   function constr = constraints(x, myauxdata)
      % Append constraints for xmin <= x <= xmax.
      xmin = myauxdata.xmin;
      xmax = myauxdata.xmax;
      
      constr = [x - xmin; 
                xmax - x(1:myauxdata.lenx_no_s)];

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

   function jvec = jacobvec(x, d, myauxdata)
      lenx_no_s = myauxdata.lenx_no_s;
      
      % Append Jacobian for x - xmin >=0.
      J = sparse(eye(length(x)));
      
      Jmax = sparse(zeros(lenx_no_s, length(x)));
      Jmax(1:lenx_no_s, 1:lenx_no_s) = -eye(lenx_no_s);

      % Append Jacobian for xmax - x (with no slacks) >= 0.
      J = [J; Jmax];
     
      dType = 0;
      if (size(d, 1) == size(x, 1))  % case: d == dx
         jvec = J * d;
         dType = 'dx';
      elseif (size(d, 1) == myauxdata.NINEQ) % case: d == dz
         dType = 'dz';
         jvec = J' * d;
      end
      
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