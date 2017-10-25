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
% g(x) = [x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2]
%
function self = MyEq()

    % y=g(x) 
    self.eval = @(x) sum(x.^2);

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