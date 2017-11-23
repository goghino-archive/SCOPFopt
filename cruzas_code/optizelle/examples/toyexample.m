% We consider:
% f(x) = @(x) sin(x(1)) + sin(x(2))
% 
% subject to
%
% g(x) = [x(1)^2 = 1;
%         x(2)^2 = 1]   
%
% h(x) = [cos(x(1)) <= 1;
%         cos(x(2)) <= 1]
%
% xmin .<= x .<= xmax
%
% Where g(x) = 0 
%       and
%       h(x) >= 0

% Optimize a simple optimization problem.
function toyexample()
   clc;
   clear all;
   close all;

    % Execute the optimization
    main();
end

% Define a simple objective.
% f(x) = sin(x(1)) + sin(x(2))
function self = MyObj(myauxdata)
   
    % Evaluation 
    self.eval = @(x) sin(x(1)) + sin(x(2));

    % Gradient
    self.grad = @(x) [cos(x(1));
                      cos(x(2))];

    % Hessian-vector product
    self.hessvec = @(x,dx) MyHessvec(x,dx, myauxdata);
                                
    %% Helper functions
   function hessvec = MyHessvec(x, dx, myauxdata)
      functionIdentifier = '[MyObj] MyHessvec';
      
      % Define Hessian matrix.
      H = [-sin(x(1)),        0; 
                    0, -sin(x(2))];
                 
      hessvec = H * dx;
                          
      if myauxdata.verbose
         %% Test for Inf, -Inf, and NaN in dx.
         if find(dx == Inf)
            fprintf('%s: Inf found in dx\n', functionIdentifier);
         end
         
         if find(dx == -Inf)
            fprintf('%s: -Inf found in dx\n', functionIdentifier);
         end
         
         if find(isnan(dx))
            fprintf('%s: NaN found in dx\n', functionIdentifier);
         end
         
         %% Test for Inf, -Inf, and NaN in jvec.
         if find(hessvec == Inf, 1)
            fprintf('%s: Inf found in hessvec\n', functionIdentifier);
         end
         
         if find(hessvec == -Inf, 1)
            fprintf('%s: -Inf found in hessvec\n', functionIdentifier);
         end
         
         if find(isnan(hessvec), 1)
            fprintf('%s: NaN found in hessvec\n', functionIdentifier);
         end
      end
   end
end

% Define a simple equality
%
% g(x) = [x(1)^2 = 1;
%         x(2)^2 = 1]         
% Optizelle  requires g(x) = 0, so we transform g(x) accordingly.
function self = MyEq(myauxdata)

    % y=g(x) 
    self.eval = @(x) [x(1)^2 - 1;
                      x(2)^2 - 1]; 

    % y=g'(x)dx (Jacobian * dx)
    self.p = @(x,dx) MyJacobvec(x, dx, myauxdata);

    % xhat=g'(x)*dy (Jacobian' * dy)
    self.ps = @(x, dy) MyJacobvec(x, dy, myauxdata);

    % xhat=(g''(x)dx)*dy (Partial derivatives of Jacobian * dx)
    self.pps = @(x,dx,dy) MyHessvec(x, dx, dy, myauxdata);
                         
    %% Helper functions.
   function jvec = MyJacobvec(x, d, myauxdata)
      functionIdentifier = '[MyEq] MyJacobvec';
      
      % Define Jacobian matrix.
      J = [2*x(1),    0;
                0, 2*x(2)];
      
      dType = 0;
      if (size(d, 1) == size(x, 1))  % case: d == dx
         dType = 'dx';
         jvec = J * d;
      elseif (size(d, 1) == myauxdata.NEQ) % case: d == dy
         jvec = J' * d;
         dType = 'dy';
      end
      
      if myauxdata.verbose
         %% Test for Inf, -Inf, and NaN in d.
         if find(d == Inf)
            fprintf('%s: Inf found in %s\n', functionIdentifier, dType);
         end
         
         if find(d == -Inf)
            fprintf('%s: -Inf found in %s\n', functionIdentifier, dType);
         end
         
         if find(isnan(d))
            fprintf('%s: NaN found in %s\n', functionIdentifier, dType);
         end
         
         %% Test for Inf, -Inf, and NaN in jvec.
         if find(jvec == Inf, 1)
            fprintf('%s: Inf found in jvec\n', functionIdentifier)
         end
         
         if find(jvec == -Inf, 1)
            fprintf('%s: -Inf found in jvec\n', functionIdentifier)
         end
         
         if find(isnan(jvec), 1)
            fprintf('%s: NaN found in jvec\n', functionIdentifier)
         end
      end
   end

   function hessvec = MyHessvec(x, dx, dy, myauxdata)
      functionIdentifier = '[MyEq] MyHessvec';
      
      % Define Hessian matrix.
      H = [2, 0; 
           0, 2];
        
      hessvec = (H * dx)' * dy; 
      
      if myauxdata.verbose
         %% Test for Inf, -Inf, and NaN in dx.
         if find(dx == Inf, 1)
            fprintf('%s: Inf found in dx\n', functionIdentifier)
         end
         
         if find(dx == -Inf, 1)
            fprintf('%s: -Inf found in dx\n', functionIdentifier)
         end
         
         if find(isnan(dx), 1)
            fprintf('%s: NaN found in dx\n', functionIdentifier)
         end
         
         %% Test for Inf, -Inf, and NaN in dy.
         if find(dy == Inf, 1)
            fprintf('%s: Inf found in dy\n', functionIdentifier)
         end
         
         if find(dy == -Inf, 1)
            fprintf('%s: -Inf found in dy\n', functionIdentifier)
         end
         
         if find(isnan(dy), 1)
            fprintf('%s: NaN found in dy\n', functionIdentifier)
         end
         
         %% Test for Inf, -Inf, and NaN in hessvec.
         if find(hessvec == Inf, 1)
            fprintf('%s: Inf found in hessvec\n', functionIdentifier)
         end
         
         if find(hessvec == -Inf, 1)
            fprintf('%s: -Inf found in hessvec\n', functionIdentifier)
         end
         
         if find(isnan(hessvec), 1)
            fprintf('%s: NaN found in hessvec\n', functionIdentifier)
         end
      end
   end
end

% Define inequalities, and bounds on x.
%
% h(x) = [cos(x(1)) <= 1;
%         cos(x(2)) <= 1;
%         x(1)      >= xmin(1);
%         x(2)      >= xmin(2);
%         x(1)      <= xmax(1);
%         x(2)      <= xmax(2);
% Optizelle  requires h(x) >= 0, so we transform g(x) accordingly.
function self = MyIneq(myauxdata)

    % z=h(x) 
    self.eval = @(x) [cos(x(1));
                      cos(x(2));
                      x - myauxdata.xmin;
                      myauxdata.xmax - x];

    % z=h'(x)dx
    self.p = @(x,dx) MyJacobvec(x, dx, myauxdata);

    % xhat=h'(x)*dz
    self.ps = @(x,dz) MyJacobvec(x, dz, myauxdata);

    % xhat=(h''(x)dx)*dz
    % Since all constraints are affine, we have h''(x) = 0.
    self.pps = @(x,dx,dz) sparse(length(x), length(x)); 
    
    %% Helper functions.
    function jvec = MyJacobvec(x, d, myauxdata)
      functionIdentifier = '[MyIneq] MyJacobvec';
      
      % Define Jacobian matrix.
      J = [  -sin(x(1)),          0;
                      0, -sin(x(2));
             sparse(eye(length(x)));
            -sparse(eye(length(x)))];
      
      dType = 0;
      if (size(d, 1) == size(x, 1))  % case: d == dx
         dType = 'dx';
         jvec = J * d;
      elseif (size(d, 1) == myauxdata.NINEQ) % case: d == dz
         jvec = J' * d;
         dType = 'd';
      end
      
      if myauxdata.verbose
         %% Test for Inf, -Inf, and NaN in d.
         if find(d == Inf)
            fprintf('%s: Inf found in %s\n', functionIdentifier, dType);
         end
         
         if find(d == -Inf)
            fprintf('%s: -Inf found in %s\n', functionIdentifier, dType);
         end
         
         if find(isnan(d))
            fprintf('%s: NaN found in %s\n', functionIdentifier, dType);
         end
         
         %% Test for Inf, -Inf, and NaN in jvec.
         if find(jvec == Inf, 1)
            fprintf('%s: Inf found in jvec\n', functionIdentifier)
         end
         
         if find(jvec == -Inf, 1)
            fprintf('%s: -Inf found in jvec\n', functionIdentifier)
         end
         
         if find(isnan(jvec), 1)
            fprintf('%s: NaN found in jvec\n', functionIdentifier)
         end
      end
   end
end

% Actually runs the program
function main()

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();
    
    % Settings...
    infiniteBounds = 0;
    bigButNotInfiniteBounds = 0;
    verbose = 1;

    NEQ = 2;   % Number of equality constraints.
    NINEQ = 6; % Number of inequality constraints.

    % Generate an initial guess 
    x0 = [-0.8; -0.8];
    
    % Define lower and upper bounds xmin .<= x .<= xmax
    if infiniteBounds
      xmin = -inf(size(x0));
      xmax = inf(size(x0));
    elseif bigButNotInfiniteBounds 
       xmin = -1e20 * ones(size(x0));
       xmax = 1e20 * ones(size(x0));
    else
      xmin = -10 * ones(size(x0));
      xmax =  10 * ones(size(x0));
    end
    
    % Allocate memory for the equality multiplier 
    y = zeros(NEQ, 1);

    % Allocate memory for the inequality multiplier 
    z = zeros(NINEQ, 1);
    
    myauxdata.xmin = xmin;
    myauxdata.xmax = xmax;
    myauxdata.NEQ = NEQ;
    myauxdata.NINEQ = NINEQ;
    myauxdata.verbose = verbose; % print out all information for testing
    
    % Create an optimization state
    state = Optizelle.Constrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,x0,y,z);

    % Create a bundle of functions
    fns = Optizelle.Constrained.Functions.t;
    fns.f = MyObj(myauxdata);
    fns.g = MyEq(myauxdata);
    fns.h = MyIneq(myauxdata);

    % Solve the optimization problem
%     tic
    state = Optizelle.Constrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout, ...
        fns,state);
%     toc
     
    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e,%e,%e)\n', state.x(1), state.x(2));
end