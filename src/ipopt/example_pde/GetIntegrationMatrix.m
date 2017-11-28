function v = GetIntegrationMatrix(IdxCube,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns a row vector to integrate a function over the whole simulation region
    % discretized according to IdxCube and h,
    % Uses center point rule for integration
    % If Y is a column vector discetization of a function y according to
    % IdxCube and h than Integral_{Omega}{y(x)dx} can be approximated by
    % Int = GetIntegrationMatrix(IdxCube,h)*y;
    % IdxCube: (tensor of order of problem dimension,
    %           size corresponding to problem discretization)
    %   e.g. for a simulation region (0,1)^3 with h1 = 0.1, h2=h3=0.2
    %   IdxCube must be a 10*5*5-tensor
    % h: (vector of length of problem dimension or single value)
    %   discretization step
    % rv: modified difference operator
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14    
    Dim = length(size(IdxCube));
    if length(h)==1 && Dim>1
        h = repmat(h,1,Dim);
    end
    v = prod(h)*ones(1,numel(IdxCube));
    str = zeros(1,4*Dim-1);
    str(:) = ' ';
    str = char(str);
    str(4*(1:Dim-1)) = ',';
    str(4*(1:Dim)-2) = ':';
    for iDim = 1:Dim
        str(4*(iDim-1)+(1:3)) = ' 1 ';
        Idx = eval(['IdxCube(' str ')']);
        v(Idx(:)) = v(Idx(:)) * 0.5;
        str(4*(iDim-1)+(1:3)) = 'end';
        Idx = eval(['IdxCube(' str ')']);
        v(Idx(:)) = v(Idx(:)) * 0.5;
        str(4*(iDim-1)+(1:3)) = ' : ';
    end
end

