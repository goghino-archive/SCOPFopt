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
function self = MyObj()
   
    % Evaluation 
    self.eval = @(x) sin(x(1)) + sin(x(2));

    % Gradient
    self.grad = @(x) [cos(x(1));
                      cos(x(2))];

    % Hessian-vector product
    self.hessvec = @(x,dx) MyHessvec(x,dx);
                                
    %% Helper functions
   function hessvec = MyHessvec(x, dx)
      functionIdentifier = '[MyObj] MyHessvec';
      
      H = [-sin(x(1)),        0; 
                    0, -sin(x(2))];
                 
      hessvec = H * dx;
                                
      %% Test for Inf, -Inf, and NaN in dx.
      if find(d == Inf)
         fprintf('%s: Inf found in dx', functionIdentifier);
      end
      
      if find(d == -Inf)
         fprintf('%s: -Inf found in dx', functionIdentifier);
      end
      
      if find(isnan(d))
         fprintf('%s: NaN found in dx', functionIdentifier);
      end
      
      %% Test for Inf, -Inf, and NaN in jvec.
      if find(hessvec == Inf, 1)
         fprintf('%s: Inf found in hessvec', functionIdentifier);
      end
      
      if find(hessvec == -Inf, 1)
         fprintf('%s: -Inf found in hessvec', functionIdentifier);
      end
      
      if find(isnan(hessvec), 1)
         fprintf('%s: NaN found in hessvec', functionIdentifier);
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
                         
    % Helper functions.
   function jacobvec = MyJacobvec(x, d, myauxdata)
      functionIdentifier = '[MyEq] MyJacobvec';
      
      % Define Jacobian matrix.
      J = [2*x(1),    0;
                0, 2*x(2)];
      
      dType = 0;
      if (size(d, 1) == size(x, 1))  % case: d == dx
         dType = 'dx';
         jacobvec = J * d;
      elseif (size(d, 1) == myauxdata.NEQ) % case: d == dy
         jacobvec = J' * d;
         dType = 'dy';
      end
      
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
         fprintf('%s: Inf found in jvec', functionIdentifier)
      end
      
      if find(jvec == -Inf, 1)
         fprintf('%s: -Inf found in jvec', functionIdentifier)
      end
      
      if find(isnan(jvec), 1)
         fprintf('%s: NaN found in jvec', functionIdentifier)
      end
   end

   function hessvec = MyHessvec(x, dx, dy, myauxdata)
      hessvec = ([2, 0; 
                  0, 2] * dx)' * dy; 
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
    self.p = @(x,dx) [  -sin(x(1)),          0;
                                 0, -sin(x(2));
                                 sparse(eye(length(x)));
                                -sparse(eye(length(x)))] * dx;

    % xhat=h'(x)*dz
    self.ps = @(x,dz) [-sin(x(1)),          0;
                                0, -sin(x(2));
                                0,          0;
                                0,          0;
                                0,          0;
                                0,          0]' * dz;

    % xhat=(h''(x)dx)*dz
    % Since all constraints are affine, we have h''(x) = 0.
    self.pps = @(x,dx,dz) sparse(length(x), length(x)); 
end

% Actually runs the program
function main()

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();
    
    % Settings.
    % Toggle on/off for xmin = -inf(size(x0)); xmax = inf(size(x0));
    infiniteBounds = 0;
    bigButNotInfiniteBounds = 1;
    
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
    
    % Create an optimization state
    state = Optizelle.Constrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,x0,y,z);

    % Create a bundle of functions
    fns = Optizelle.Constrained.Functions.t;
    fns.f = MyObj();
    fns.g = MyEq(myauxdata);
    fns.h = MyIneq(myauxdata);

    % Solve the optimization problem
    tic;
    state = Optizelle.Constrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout, ...
        fns,state);
    myTime = toc
     
    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e,%e,%e)\n', state.x(1), state.x(2));
end