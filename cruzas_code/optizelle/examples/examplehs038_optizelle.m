% Optimize a simple optimization problem.
function examplehs038_optizelle()
    % Execute the optimization
    main();
end

% Define a simple objective.
function self = MyObj()

    % Evaluation 
    self.eval = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);

    % Gradient
    self.grad = @(x) [
       -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1)); 
       200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1); 
       -360*x(3)*(x(4)-x(3)^2) - 2*(1-x(3)); 
       180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)];

    % Hessian-vector product
    self.hessvec = @(x,dx) [
       (1200*x(1)^2-400*x(2)+2)*dx(1);
        (-400*x(1))*dx(1) + 220.2*dx(2);
        (1080*x(3)^2 - 360*x(4) + 2)*dx(3);
        19.8*dx(2) - (360*x(3))*dx(3) + 200.2*dx(4)];
end
   
% Actually runs the program
function main()

    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Generate an initial guess 
%     x = [-3, -1, -3, -1]; % This is not working...
    x = [-3; -1; -3; -1];

    % Create an optimization state
    state = Optizelle.Unconstrained.State.t( ...
        Optizelle.Rm, x);

     fname = '/Users/samuelcruz/Documents/GitHub/SCOPFopt/cruzas_code/tr_newton.json';
    % Read the parameters from file
    state = Optizelle.json.Unconstrained.read( ...
        Optizelle.Rm,fname,state);

    % Create a bundle of functions
    fns = Optizelle.Unconstrained.Functions.t;
    fns.f = MyObj();

    % Solve the optimization problem
    state = Optizelle.Unconstrained.Algorithms.getMin( ...
        Optizelle.Rm, Optizelle.Messaging.stdout, ...
        fns,state);

    % Print out the reason for convergence
    fprintf('The algorithm converged due to: %s\n', ...
        Optizelle.OptimizationStop.to_string(state.opt_stop));

    % Print out the final answer
    fprintf('The optimal point is: (%e,%e,%e,%e)\n', state.x(1), state.x(2), ...
                                                     state.x(3), state.x(4));

    % Write out the final answer to file
%     Optizelle.json.Constrained.write_restart(Optizelle.Rm, 'solution.json', ...
%                                              state);
end