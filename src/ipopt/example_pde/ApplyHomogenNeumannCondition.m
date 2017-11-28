function Op = ApplyHomogenNeumannCondition(Op,IdxCube,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Applies homogenious Neumann boundary condition (BC)
	% dy/dn = 0
    % at any border point in IdxCube to the difference operator Op,
    % uses cerntral difference quotient to approximate dy/dn
    % Op: (Difference Operator e.g. as received from GetStandardOp)
    % IdxCube: (tensor of order of problem dimension,
    %           size corresponding to problem discretization)
    %   e.g. for a simulation region (0,1)^3 with h1 = 0.1, h2=h3=0.2
    %   IdxCube must be a 10*5*5-tensor
    % h: (vector of length of problem dimension or single value)
    %   discretization step
    % return value:
    % Op: modified difference operator
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14

    Dim = length(size(IdxCube));
    if length(h) == 1 && Dim>1
        h = repmat(h,1,Dim);
    end

    str = char(' ' * ones(1,4*Dim-1));
    str(4*(1:Dim-1)) = ',';
    for iDim=1:Dim
    	str(4*(1:Dim)-2) = '1';
        str((4*(iDim-1))+(1:3)) = ' 2 ';
        IdxStep = eval(['IdxCube(' str ')']) - 1;
	
        str(4*(1:Dim)-2) = ':';
        str((4*(iDim-1))+(1:3)) = ' 1 ';
        Idx = eval(['IdxCube(' str ')']);
        Op = ApplyBoundaryCondition(Idx,-IdxStep,1,0,0,Op,0,h(iDim));

        str((4*(iDim-1))+(1:3)) = 'end';
        Idx = eval(['IdxCube(' str ')']);
        Op = ApplyBoundaryCondition(Idx, IdxStep,1,0,0,Op,0,h(iDim));

        str((4*(iDim-1))+(1:3)) = ' : ';
    end
end
