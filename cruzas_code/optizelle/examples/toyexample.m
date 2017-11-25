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
    self.grad = @(x) MyGrad(x, myauxdata);

    % Hessian-vector product
    self.hessvec = @(x,dx) MyHessvec(x,dx, myauxdata);
                                
    %% Helper functions
   function grad = MyGrad(x, myauxdata)
      grad = [cos(x(1));
              cos(x(2))];
      
      if myauxdata.withSlacks
         grad = [grad; 0; 0];
      end
   end
    
   function hessvec = MyHessvec(x, dx, myauxdata)
      functionIdentifier = '[MyObj] MyHessvec';
      
      % Define Hessian matrix.
      H = [-sin(x(1)),        0; 
                    0, -sin(x(2))];
                 
      if myauxdata.withSlacks
         tmp_H = sparse(zeros(length(x)));
         
         tmp_H(1:2, 1:2) = H;
         
         H = tmp_H;
      end
                 
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
% 
% If using slack variables, we have:
% g(x) = [x(1)^2 = 1;
%         x(2)^2 = 1;
%         cos(x(1)) + s = 1;
%         cos(x(2)) + s = 1;] 
function self = MyEq(myauxdata)

    % y=g(x) 
    self.eval = @(x) MyEval(x, myauxdata);

    % y=g'(x)dx (Jacobian * dx)
    self.p = @(x,dx) MyJacobvec(x, dx, myauxdata);

    % xhat=g'(x)*dy (Jacobian' * dy)
    self.ps = @(x, dy) MyJacobvec(x, dy, myauxdata);

    % xhat=(g''(x)dx)*dy (Partial derivatives of Jacobian * dx)
    self.pps = @(x,dx,dy) MyHessvec(x, dx, dy, myauxdata);
                         
    %% Helper functions.
   function res = MyEval(x, myauxdata)
      gx = [x(1)^2 - 1;
             x(2)^2 - 1];
        
      if myauxdata.withSlacks
         hx = [1 - cos(x(1));
               1 - cos(x(2))];
            
         % Extract slack variables.
         s = [x(3); x(4)];
            
         res = [gx; hx - s];
      else
         res = gx;
      end
   end
    
   function jvec = MyJacobvec(x, d, myauxdata)
      functionIdentifier = '[MyEq] MyJacobvec';
      
      % Define Jacobian matrix.
      J = [2*x(1),      0;
                0, 2*x(2)];
             
      if myauxdata.withSlacks 
         tmp_J = sparse(zeros(length(x)));
         
         tmp_J(1:2, 1:2) = J;
         
         J = tmp_J;
         
         dhx = [sin(x(1)),          0,  -1,  0;
                      0,    sin(x(2)),  0,  -1];
         J = [J; dhx];
      end
      
      dType = 0;
      if (size(d, 1) == size(x, 1))  % case: d == dx
         dType = 'dx';
         disp(dType)

         jvec = J * d;
      elseif (size(d, 1) == myauxdata.NEQ) % case: d == dy
         dType = 'dy';
         disp(dType)
         
         jvec = J' * d;
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
        
      if myauxdata.withSlacks
         tmp_H = sparse(zeros(length(x)));
         
         ddhx = [cos(x(1)),      0;
                        0, cos(x(2))];
         
         tmp_H(1:2, 1:2) = H + ddhx;
         
         H = tmp_H;
      end
        
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
% Note, the "inequality constraints" xmin .<= x .<= xmax.
function self = MyIneq(myauxdata)

    % z=h(x) 
    self.eval = @(x) MyEval(x, myauxdata);

    % z=h'(x)dx
    self.p = @(x,dx) MyJacobvec(x, dx, myauxdata);

    % xhat=h'(x)*dz
    self.ps = @(x,dz) MyJacobvec(x, dz, myauxdata);

    % xhat=(h''(x)dx)*dz
    % Since all constraints are affine, we have h''(x) = 0.
    self.pps = @(x,dx,dz) sparse(length(x), length(x)); 
    
    %% Helper functions.
   function res = MyEval(x, myauxdata)
      
      if myauxdata.withSlacks
         res = [x - myauxdata.xmin;
                myauxdata.xmax - [x(1); x(2)]];
      else
            res = [1 - cos(x(1));
                   1 - cos(x(2));
                   x - myauxdata.xmin;
                   myauxdata.xmax - x];
      end
   end
    
    function jvec = MyJacobvec(x, d, myauxdata)
      functionIdentifier = '[MyIneq] MyJacobvec';
      
      if myauxdata.withSlacks
         % Jacobian for x - xmin.
         Jmin = sparse(eye(length(x)));
         
         % Jacobian for xmax - x_no_s.
         Jmax = sparse(zeros(2, length(x)));
         
         Jmax(1:2, 1:2) = -eye(2);
         
         J = [Jmin; Jmax];
              
      else
         % Define Jacobian matrix.
         J = [  sin(x(1)),          0;
                        0, sin(x(2));
                sparse(eye(length(x)));
               -sparse(eye(length(x)))];
      end
    
      dType = 0;
      if (size(d, 1) == size(x, 1))  % case: d == dx
         dType = 'dx';
         jvec = J * d;
      elseif (size(d, 1) == myauxdata.NINEQ) % case: d == dz
         jvec = J' * d;
         dType = 'dz';
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
    % infiniteBounds and bigButNotInfiniteBounds can be left as 0, as
    % bounds of [-10; -10] for xmin and [10; 10] for xmax will be chosen in
    % that case.
    infiniteBounds = 0; 
    bigButNotInfiniteBounds = 1; 
    withSlacks = 1;
    verbose = 0;  % print message when Inf, -Inf, or NaN found
    myTol = 1e-12; % tolerance value for assertion of solution

    if infiniteBounds
      disp('Using infinite bounds...')
    elseif bigButNotInfiniteBounds
       disp('Using big but not infinite bounds...')
    end
   
    % Generate an initial guess.
    x0 = [-0.4; -0.4];
   
    % Slack variables, corresponding to number of inequality constraints.
    s = zeros(2, 1);
    
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
    
    if withSlacks
       % Note, no upper bounds on slack variables, as per Optizelle
       % documentation.
       x = [x0; s];
       xmin = [xmin; zeros(size(s))];
    else
       x = x0;
    end
    
    if withSlacks
       NEQ = 4
       NINEQ = 2*length(x0) + length(s)
    else
      NEQ = 2   % Number of equality constraints.
      NINEQ = 2*length(x0) + 2 % Number of inequality constraints.
    end
    
    % Allocate memory for the equality multiplier 
    y = zeros(NEQ, 1);
    % Allocate memory for the inequality multiplier 
    z = zeros(NINEQ, 1);
    
    x_size = size(x)
    y_size = size(y)
    z_size = size(z)
    
    myauxdata.xmin = xmin;
    myauxdata.xmax = xmax;
    myauxdata.NEQ = NEQ;
    myauxdata.NINEQ = NINEQ;
    myauxdata.withSlacks = withSlacks;
    myauxdata.verbose = verbose; % print out all information for testing
    
    % Create an optimization state
    state = Optizelle.Constrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,x,y,z);

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
    
    % Assert that all equality constraints are satisfied
    assert(abs(state.x(1)^2 - 1) < myTol)
    assert(abs(state.x(2)^2 - 1) < myTol)
    
    % Assert that all inequality constraints are satisfied
    assert(cos(state.x(1)) <= 1);
    assert(cos(state.x(2)) <= 1);
    assert(state.x(1) >= xmin(1))
    assert(state.x(2) >= xmin(2));
    assert(state.x(1) <= xmax(1))
    assert(state.x(2) <= xmax(2));
end