function examplehs051_optizelle()
    % Execute the optimization
    main();
end

% Define a simple objective.
function self = MyObj()

    % Evaluation 
    self.eval = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + ...
                     (x(5) - 1)^2;

    % Gradient
    self.grad = @(x) 2*[ x(1) - x(2);
	  x(2) + x(3) - 2 - x(1) + x(2);
	  x(2) + x(3) - 2;
	  x(4) - 1;
	  x(5) - 1 ];

    % Hessian-vector product
    self.hessvec = @(x,dx) hessvec(x, dx);
end

function H_dx = hessvec(x, dx)
    % Allocate memory for the dense Hessian in packed storage
    H = zeros(5);
    
    H(1, :) = [2, -2, 0, 0, 0];
    H(2, :) = [ -2,  4, 2, 0, 0];
    H(3, :) = [0,  2, 2, 0, 0];
    H(4, :) = [0,  0, 0, 2, 0];
    H(5, :) = [0,  0, 0, 0, 2];

    % Compute the Hessian-vector product
    H_dx = zeros(5, 1);
    for i = 1:5 
        for j = 1:5 
            H_dx(i) = H_dx(i) + H(i, j)*dx(j);
        end
    end
end

%
% g(x)= [ x(1) + 3*x(2) ]
%       [ x2 x3 - 5 x4 x5                       ]
%       [ x1^3 + x2^3 + 1                       ]
%
function self = MyEq()

    % y=g(x) 
    self.eval = @(x) [
        (sq(x(itok(1))) + sq(x(itok(2))) + sq(x(itok(3))) ...
            + sq(x(itok(4))) + sq(x(itok(5))) - 10.);
        (x(itok(2))*x(itok(3)) - 5.*x(itok(4))*x(itok(5)));
        (pow(x(itok(1)),3) + pow(x(itok(2)),3) + 1.)];

    % y=g'(x)dx
    self.p = @(x,dx) reshape(generateJac(x),3,5)*dx; 

    % xhat=g'(x)*dy
    self.ps = @(x,dy) reshape(generateJac(x),3,5)'*dy; 

    % xhat=(g''(x)dx)*dy
    self.pps = @(x,dx,dy) pps(x,dx,dy); 
end
  