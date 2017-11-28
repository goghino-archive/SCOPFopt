function W = GetInterpolationMatrix(P,varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns a matrix with linera interpolation weights at points declared
    % in P
    % P: m*n matrix for m point each with n coordinates
    % vargin: variable number of arguments according to the definition
    %   region dimension as returned by ndgrid
    %       e.g. in 3D: P = GetMeasureCoords;
    %                   [x1,x2,x3] = ndgrid(0:h1:1,0:h2:2,0:h3:0.5);
    %                   W = GetInterpolationMatrix(P,x1,x2,x3);
    % return value:
    % W: m*k - Matrix, where k is the number of discretization points
    %   Linear interpolation of a discretized funktion y (column vector)
    %   is than obtained by W*y
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14    
    Dim = size(P,2);
    
    if length(varargin)~=Dim
        error('wrong number of argument')
    end
    
    m = size(P,1);
    s = numel(varargin{1});
    h = zeros(1,Dim);
    UsedSmplPts = true( [m s] );
    dx = cell(size(varargin));
    
    str1 = char(' ' * ones(1,2*Dim-1));
    str1(2*(1:Dim-1)) = ',';
    str1(2*(1:Dim)-1) = '1';
    str2 = str1;
    str2(2*(1:Dim)-1) = '2';
    
    for iDim=1:Dim
        dx{iDim} = abs(repmat(P(:,iDim),1,s) - repmat(reshape(varargin{iDim},1,[]),m,1));
        h(iDim) = eval(['varargin{iDim}(' str2 ') - varargin{iDim}(' str1 ');']);
        UsedSmplPts = UsedSmplPts & dx{iDim}<h(iDim);
    end
    W = sparse(m,s);
    W(UsedSmplPts) = 1;
    for iDim=1:Dim
        W(UsedSmplPts) = W(UsedSmplPts).*(1-dx{iDim}(UsedSmplPts)/h(iDim));
    end
end
