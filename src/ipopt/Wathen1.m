function [status, y, u, q, IPOPT_Statistic] = Wathen1(Dim,NumPtsPerDim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Wathen1
%   Example program for usage of the Ipopt package in matlab
%   We solve the PDE constrained problem 1 in Rees/Dollar/Wathen:
%   "Optimal solvers for PDE-constrained optimization",
%   RAL-TR-2008-018, ISSN 1358-6254
%   The problem is:
%   Minimize f(y,u) = int_{Omega}{1/2 (y-y_d)^2 + alpha u^2 dx}
%   subject to  -Laplacian(y) = u          		in Omega
%               y=y_d                               on Border(Omega)
%   where Omage = (0,1)^Dim, Dim + 2,3 and
%	y_d = (2*x1-1).^2*(2*x2-1).^2			in (0,0.5)^2, 0 otherwise for Dim = 2
%	y_d = (2*x1-1).^2*(2*x2-1).^2*(2*x3-1).^2	in (0,0.5)^3, 0 otherwise for Dim = 3
%   alpha = 0.001.
%   We discretize with Finite Difference.
%   Author: Johannes Huber
%           Mathematische Institut
%           University of Basel
%   2009_01_14

% Initialize some variables to setup the problem in callback functions
% Since we work with nested functions we can precompute some frequently
% used values here and use it in the function
alpha = 0.01;
h = 1/(NumPtsPerDim-1);                             % grid step size in all directions
MeshData = ones(Dim,3);
MeshData(:,1) = 0;
MeshData(:,2) = h;
switch(Dim)
   case 2
      [x1,x2,IdxCube] = CreateMesh(MeshData);
      yd = CalcForcingFunc(x1,x2);
      IdxInner = IdxCube(2:end-1,2:end-1);
   case 3
      [x1,x2,x3,IdxCube] = CreateMesh(MeshData);
      yd = CalcForcingFunc(x1,x2,x3);
      IdxInner = IdxCube(2:end-1,2:end-1,2:end-1);
   otherwise
      error('Dimension not implemented');
end
yd = yd(:);
NumPts = numel(IdxCube);
% Compute the indices in inner grid points and at the border
IdxInner = IdxInner(:);
IdxBound = setdiff(IdxCube(:),IdxInner);
% setup for ipopt, i.e. initial values and bounds
y0 = yd';                                            % starting guess for discretized function y
u0 = ones(1,NumPts);                                % starting guess for discretized function u
v0 = [y0 u0];                                       % Unknown variables are y and u, combined as v
lb = [-Inf *ones(size(y0)) -Inf*ones(size(u0))];     % Lower bound of variables
ub = [ Inf *ones(size(y0))  Inf*ones(size(u0))];     % Upper bound of variables
% In every grid point we must fulfill the discretized PDE which are the
% equality constraints.
lbc = zeros(1,NumPts);
ubc = zeros(1,NumPts);

IntM = GetIntegrationMatrix(IdxCube,h);
SysOp = GetStandardOp('Laplace',IdxCube,h);
SysOp = ApplyStencil(IdxCube,-1,IdxBound,SysOp);
SysM = -CreateSystemMatrix(SysOp);

OutputFilename = ['Wathen_1_D' num2str(Dim) '_' num2str(NumPtsPerDim) '.txt'];

options.lb = lb;  % Lower bound on the variables.
options.ub = ub;  % Upper bound on the variables.
options.cl = lbc;   % Lower bounds on the constraint functions.
options.cu = ubc;   % Upper bounds on the constraint functions.

% The callback functions.
funcs.objective         = @computeObjective;
funcs.constraints       = @computeConstraints;
funcs.gradient          = @computeGradient;
funcs.jacobian          = @computeJacobian;
funcs.jacobianstructure = @;
funcs.hessian           = @computeHessian;
funcs.hessianstructure  = @;

if NumPtsPerDim<5
   % It's posible to ask ipopt to doublecheck the derivatives.
   % for small problem sizes we do that
   [status v multipliers] = ipopt(v0,lb,ub,lbc,ubc,@computeObjective,...
      @computeGradient,@computeConstraints,...
      @computeJacobian,@computeHessian,[],'',[],...
      'derivative_test','second-order',...
      'max_iter',10,...
      'output_file',OutputFilename);
else
   [status v multipliers] = ipopt(v0,lb,ub,lbc,ubc,@computeObjective,...
      @computeGradient,@computeConstraints,...
      @computeJacobian,@computeHessian,[],'',[],...
      'obj_scaling_factor',1/h^2,...
      'max_iter',500,...
      'output_file',OutputFilename);
end

% Get IPOPT output data
IPOPT_Statistic = ReadIPOPTOutputFile(OutputFilename);
delete(OutputFilename);

v = reshape(v,NumPts,[]);
switch(Dim)
   case 2
      y = reshape(v(:,1),NumPtsPerDim,NumPtsPerDim);
      u = reshape(v(:,2),NumPtsPerDim,NumPtsPerDim);
      q = reshape(multipliers.lambda,NumPtsPerDim,NumPtsPerDim);
   case 3
      y = reshape(v(:,1),NumPtsPerDim,NumPtsPerDim,NumPtsPerDim);
      u = reshape(v(:,2),NumPtsPerDim,NumPtsPerDim,NumPtsPerDim);
      q = reshape(multipliers.lambda,NumPtsPerDim,NumPtsPerDim,NumPtsPerDim);
   otherwise
      error('Dimension not implemented');
end
%   End of main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Start of nested callback functions
%   ipopt will request these functions for function evaluation
%   recieve current y and u values in v=[y u] as row vector
   function f = computeObjective(v)
      % Computes objective function i.e. function that has to be minimized
      v = reshape(v,NumPts,[]);
      tmp = 0.5*(v(:,1)-yd).^2+alpha*v(:,2).^2;
      f = IntM*tmp;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function g = computeGradient (v)
      % Computes gradient of the objective function as a column vector
      v = reshape(v,NumPts,[]);
      g = v;
      g(:,1) = IntM.*(v(:,1)-yd)';
      g(:,2) = 2*alpha*IntM.*v(:,2)';
      g = reshape(g,[],1);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function c = computeConstraints (v)
      % Computes the constraint functions as column vector
      % i.e. in every grid point the error value for the discretized PDE
      % -Laplace(y) = u	in Omega
      % y = yd		on Border(Omega)
      % Every row belongs to one equality and thus to one grid point
      v = reshape(v,NumPts,[]);
      c = SysM*v(:,1);
      c(IdxInner) = c(IdxInner)-v(IdxInner,2);
      c(IdxBound) = c(IdxBound)-yd(IdxBound);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function J = computeJacobian (v, returnStructureOnly)
      % Computes the Jacobian of the constraint function
      if returnStructureOnly
         % ipopt asks for the sparsity pattern of the Jaconbian matrix
         % i.e. for the matrix entries that are potentially nonzero
         % Set those entries here to one, all others to zero
         J = [SysM -speye(NumPts)];
      else
         % Compute the actual Jacobian
         tmp = zeros(NumPts,1);
         tmp(IdxInner) = -1;
         J = [SysM spdiags(tmp,0,NumPts,NumPts)];
      end
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function computeHessian is not neccesary, if Hessian is
%   approxiamted by ipopt, but it is advisable to implement it if its
%   evaluation is not too expensive. Ipopt then will converge faster
%   and will be more robust.
   function H = computeHessian(v, sigma, lambda, returnStructureOnly)
      % Computes the Hessian matrix of the function
      % L = sigma*f+lamda*c
      % where f is the objective function and c is the (vector values)
      % constrained function
      % This function is similar to the Lagrangian except the factor
      % sigma which exists for masking f
      if returnStructureOnly
         % ipopt asks for the sparsity pattern of the Hessian matrix
         % i.e. for the matrix entries that are potentially nonzero
         % Set those entries here to one, all others to zero
         H = speye(2*NumPts,2*NumPts);
      else
         % Compute the actual Hessian matrix
         Diag = [IntM';2*alpha*IntM'];
         H = spdiags(sigma*Diag,0,2*NumPts,2*NumPts);
      end
   end
%   End of nested callback functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%   End of function Wathen1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rv = CalcForcingFunc(x1,x2,x3)
rv = zeros(size(x1));
if exist('x3')
   idx1 = find(x1(:,1,1)<=0.5);
   idx2 = find(x2(1,:,1)<=0.5);
   idx3 = find(x3(1,1,:)<=0.5);
   rv(idx1,idx2,idx3) = (2*x1(idx1,idx2,idx3)-1).^2.*(2*x2(idx1,idx2,idx3)-1).^2.*(2*x3(idx1,idx2,idx3)-1).^2;
else
   idx1 = find(x1(:,1)<=0.5);
   idx2 = find(x2(1,:)<=0.5);
   rv(idx1,idx2) = (2*x1(idx1,idx2)-1).^2.*(2*x2(idx1,idx2)-1).^2;
end
end
