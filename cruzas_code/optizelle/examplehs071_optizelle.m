% Optimize a simple optimization problem.
function examplehs071_optizelle()
    % Execute the optimization
    main();
end

% Define a simple objective.
function self = MyObj()

    % Evaluation 
    self.eval = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);

    % Gradient
    self.grad = @(x) [x(1)*x(4) + x(4)*sum(x(1:3));
        x(1)*x(4);
        x(1)*x(4) + 1;
        x(1)*sum(x(1:3))];

    % Hessian-vector product
    self.hessvec = @(x,dx) hessvec(x, dx);
end

function H_dx = hessvec(x, dx)
    % Allocate memory for the dense Hessian in packed storage
    H = zeros(5);
    
    H(1, :) = [2*x(4), x(4), x(4), 2*x(1) + x(2) + x(3)];
    H(2, :) = [x(4),  0,  0, x(1)];
    H(3, :) = [x(4),  0,  0, x(1)];
    H(4, :) = [2*x(1) + x(2) + x(3), x(1), x(1), 0];
    
    % Compute the Hessian-vector product
    % REVIEW: can we simply change this to ... H_dx = H*dx ?
    H_dx = zeros( size(H, 1), 1); 
    for i = 1:size(H, 1)
        for j = 1:size(H, 1)
            H_dx(i) = H_dx(i) + H(i, j)*dx(j);
        end
    end
end

% Define a simple equality
%
% g(x) = [ x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 = 40]
%
function self = MyEq()

    % y=g(x) 
    self.eval = @(x) sum(x.^2) - 40;

    % y=g'(x)dx
    self.p = @(x,dx) [2*x(1)*dx(1) + 2*x(2)*dx(2) + ...
                      2*x(3)*dx(3) + 2*x(4)*dx(4)];

    % xhat=g'(x)*dy
    self.ps = @(x,dy)  [2*x(1)*dy(1);
                        2*x(2)*dy(2);
                        2*x(3)*dy(3);
                        2*x(4)*dy(4)];

    % xhat=(g''(x)dx)*dy
    self.pps = @(x,dx,dy) [2*dx(1)*dy(1);
                           2*dx(2)*dy(2);
                           2*dx(3)*dy(3);
                           2*dx(4)*dy(4)]; 
end

% Define inequalities, and bounds on x
%
% h(x) = [ x(1)*x(2)*x(3)(x(4) >= 25 ] 
%        [ x(1) >= 1]
%        [ x(2) >= 1]
%        [ x(3) >= 1]
%        [ x(4) >= 1]
%        [x(1) <= 5] = [ -x(1) >= -5]
%        [x(2) <= 5] = [ -x(2) >= -5]
%        [x(3) <= 5] = [ -x(3) >= -5]
%        [x(4) <= 5] = [ -x(4) >= -5]
% We note that expressing x(i) >= 1, given that n=4 means that sum(x) >= 4.
% Similarly, since -x(i) >= -5, given that n=4, -sum(x) >= -20. 
% So we have:
% h(x) = [ x(1)*x(2)*x(3)(x(4) >= 25 ] 
%        [ x(1) + x(2) + x(3) + x(4) >= 4]
%        [ -(x(1) + x(2) + x(3) + x(4)) >= -20]
function self = MyIneq()

    % z=h(x) 
    self.eval = @(x) [prod(x) - 25; 
                      sum(x) - 4;
                      -sum(x) + 20];

    % z=h'(x)dx
    self.p = @(x,dx) generateJac(x)*dx;

    % xhat=h'(x)*dz
    self.ps = @(x,dz) generateJac(x)' * dy;

    % xhat=(h''(x)dx)*dz
    self.pps = @(x,dx,dz) [ 0. ]; 
end

% Generate a dense version of the Jacobian
function jac = generateJac(x)
   % Jacobian dimension is: number of constraints by number of variables.
   jac = zeros(3, 4);
   
   % First row corresponds to partial derivatives of first constraint.
   jac(1,1) = x(2)*x(3)*x(4);
   jac(1,2) = x(1)*x(3)*x(4);
   jac(1,3) = x(1)*x(2)*x(4);
   jac(1,4) = x(1)*x(2)*x(3);
   
   % Second row corresponds to partial derivatives of second constraint.
   jac(2,1) = 1;
   jac(2,2) = 1;
   jac(2,3) = 1;
   jac(2,4) = 1;
   
   % Third row corresponds to partial derivatives of third constraint.
   jac(3,1) = -1;
   jac(3,2) = -1;
   jac(3,3) = -1;
   jac(3,4) = -1;
end