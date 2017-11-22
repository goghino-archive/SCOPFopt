% Optimize a simple optimization problem.
function toyexample()
   clc;
   clear all;
   close all;

    % Execute the optimization
    main();
end

% Define a simple objective.
% f(x) = x^2 + y^2
function self = MyObj()

    % Evaluation 
    self.eval = @(x) x(1)^2 + x(2)^2;

    % Gradient
    self.grad = @(x) [2*x(1);
                      2*x(2)];

    % Hessian-vector product
    self.hessvec = @(x,dx) [2, 0; 
                            0, 2] * dx;
end

% Define a simple equality
%
% g(x) = (x(1) - 2)^2 + (x(2) - 2)^2 = 1
%
function self = MyEq()

    % y=g(x) 
    self.eval = @(x) (x(1) - 2)^2 + (x(2) - 2)^2 - 1; 

    % y=g'(x)dx
%     self.p = @(x,dx) [2*(x(1) - 2), 2*(x(2) - 2)] * dx;
   self.p = @(x,dx) bs(x, dx);

   function res = bs(x, dx)
%       x_size = size(x)
%       dx_size = size(dx)
      
      res = [2*(x(1) - 2), 2*(x(2) - 2)] * dx;
   end
    
    % xhat=g'(x)*dy
    self.ps = @(x,dy)  [2*(x(1) - 2);
                        2*x(2)] .* dy;

    % xhat=(g''(x)dx)*dy
    self.pps = @(x,dx,dy) ([2, 0; 
                            0, 2] * dx) .* dy; 
end

% Define inequalities, and bounds on x.
%
% h(x) = [sin(x(1)) <= 1; 
%         sin(x(2)) <= 1;]
function self = MyIneq()

    % z=h(x) 
    self.eval = @(x) [1 - sin(x(1));
                      1 - sin(x(2))];

    % z=h'(x)dx
    self.p = @(x,dx) [-cos(x(1)), 0;
                      0         , -cos(x(2))] * dx;

    % xhat=h'(x)*dz
    self.ps = @(x,dz) [-cos(x(1)), 0;
                      0         , -cos(x(2))] .* dz;

    % xhat=(h''(x)dx)*dz
    % Since all constraints are affine, we have h''(x) = 0.
    self.pps = @(x,dx,dz) ([sin(x(1)), 0
                           0        , sin(x(2))] * dx) .* dz; 
end

% Actually runs the program
function main()

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Generate an initial guess 
    x = [0; 0];
    
    % Allocate memory for the equality multiplier 
    y = zeros(1, 1);

    % Allocate memory for the inequality multiplier 
    z = zeros(2, 1);
    
    % Create an optimization state
    state = Optizelle.Constrained.State.t( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,x,y,z);

    % Create a bundle of functions
    fns = Optizelle.Constrained.Functions.t;
    fns.f = MyObj();
    fns.g = MyEq();
    fns.h = MyIneq();

    % Solve the optimization problem
    state = Optizelle.Constrained.Algorithms.getMin( ...
        Optizelle.Rm,Optizelle.Rm,Optizelle.Rm,Optizelle.Messaging.stdout, ...
        fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e,%e,%e)\n', state.x(1), state.x(2));
end