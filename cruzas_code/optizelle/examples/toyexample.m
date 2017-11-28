% We consider:
% f(x) = @(x) sin(x(1)) + sin(x(2))
% 
% subject to
%
% g(x) = [x(1)^2 - 1 = 0;
%         x(2)^2 - 1 = 0]   
%
% h(x) = [1 - cos(x(1)) >= 0;
%         1 - cos(x(2)) >= 0]
%
% x - xmin >= 0
% xmax - x >= 0 

% Optimize a simple optimization problem.
function toyexample(fname)
    if nargin < 1
        fname = '';
    end

    % Execute the optimization
    main(fname);
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
              cos(x(2));
              0;
              0];
   end
    
   function hessvec = MyHessvec(x, dx, myauxdata)
      functionIdentifier = '[MyObj] MyHessvec';
      
      % Define Hessian matrix.
      H = [-sin(x(1)), 0, 0, 0; 
           0, -sin(x(2)), 0, 0;
           0, 0, 0, 0;
           0, 0, 0, 0];
                 
      hessvec = H * dx;
   end
end

% Define a simple equality
%
% g(x) = [x(1)^2 - 1 = 0;
%         x(2)^2 - 1 = 0]   
%
% h(x) = [1 - cos(x(1)) - s = 0;
%         1 - cos(x(2)) - s = 0]
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
      res = [x(1)^2 - 1;
             x(2)^2 - 1;
             1 - cos(x(1)) - x(3);
             1 - cos(x(2)) - x(4)];
   end
    
   function jvec = MyJacobvec(x, d, myauxdata)
      functionIdentifier = '[MyEq] MyJacobvec';
      
      % Define Jacobian matrix.
      J = [2*x(1), 0, 0, 0;
           0, 2*x(2), 0, 0;
           sin(x(1)), 0, -1,  0;
           0, sin(x(2)), 0,  -1];
             
      if (size(d, 1) == size(x, 1))  % case: d == dx
         jvec = J * d;
      elseif (size(d, 1) == myauxdata.NEQ) % case: d == dy   
         jvec = J' * d; % jvec = J .* d;
      end
   end

   function hessvec = MyHessvec(x, dx, dy, myauxdata)
      functionIdentifier = '[MyEq] MyHessvec';
      
      % Define Hessian matrix.
      H = [2*dy(1), 0, 0, 0; 
           0, 2*dy(2), 0, 0;
           0, 0, 0, 0;
           0, 0, 0, 0] ...
       +  [cos(x(1))*dy(3), 0, 0, 0;
           0, cos(x(2))*dy(4), 0, 0; 
           0, 0, 0, 0;
           0, 0, 0, 0];
        
      hessvec = H * dx; 
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
      
     res = [x - myauxdata.xmin;
            myauxdata.xmax - [x(1); x(2)]];
   end
    
    function jvec = MyJacobvec(x, d, myauxdata)
      functionIdentifier = '[MyIneq] MyJacobvec';
         
      J = [1,0,0,0;
           0,1,0,0;
           0,0,1,0;
           0,0,0,1;
           -1,0,0,0;
           0,-1,0,0;];
    
      if (size(d, 1) == size(x, 1))  % case: d == dx
         jvec = J * d;
      elseif (size(d, 1) == myauxdata.NINEQ) % case: d == dz
         jvec = J' * d; %jvec = J .* d;
      end
   end
end

% Actually runs the program
function main(fname)

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();
    
    % Settings...
    % infiniteBounds and bigButNotInfiniteBounds can be left as 0, as
    % bounds of [-10; -10] for xmin and [10; 10] for xmax will be chosen in
    % that case.
    infiniteBounds = 0; 
    bigButNotInfiniteBounds = 0; 
    myTol = 1e-6; % tolerance value for assertion of solution

    if infiniteBounds
      disp('Using infinite bounds...')
    elseif bigButNotInfiniteBounds
       disp('Using big but not infinite bounds...')
    end
   
    % Generate an initial guess.
    x0 = [-0.5; -0.5];
   
    % Slack variables, corresponding to number of inequality constraints.
    s = zeros(2, 1) + 1e-5;
    smin = zeros(2, 1);
    
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
    
   x = [x0; s];
   xmin = [xmin; smin];
    
   NEQ = 4;
   NINEQ = 2*length(x0) + length(s);
    
    % Allocate memory for the equality multiplier 
    y = ones(NEQ, 1);
    
    % Allocate memory for the inequality multiplier 
    z = ones(NINEQ, 1);
    
    myauxdata.xmin = xmin;
    myauxdata.xmax = xmax;
    myauxdata.NEQ = NEQ;
    myauxdata.NINEQ = NINEQ;
    
    % Create an optimization state
    state = Optizelle.Constrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,x,y,z);
    
    % Read the parameters from file
    if(~strcmp(fname, ''))
    state = Optizelle.json.Constrained.read( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,fname,state);
    end
    
    %Modify the state so that we just run our diagnostics and exit
%     state.dscheme = Optizelle.DiagnosticScheme.DiagnosticsOnly;
%     state.f_diag = Optizelle.FunctionDiagnostics.SecondOrder;
%     state.x_diag = Optizelle.VectorSpaceDiagnostics.Basic;
%     state.h_diag = Optizelle.FunctionDiagnostics.SecondOrder;
%     state.z_diag = Optizelle.VectorSpaceDiagnostics.EuclideanJordan;
%     state.L_diag = Optizelle.FunctionDiagnostics.SecondOrder;

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
    fprintf('The optimal point is: (%e,%e,%e,%e)\n', state.x(1), state.x(2), state.x(3), state.x(4));
    
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