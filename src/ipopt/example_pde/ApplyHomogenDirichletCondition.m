function [Op, rhs] = ApplyHomogenDirichletCondition(Op,IdxCube,rhs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Applies homogenious Dirichlet boundary condition (BC)
	% y = 0
    % at any border point in IdxCube to the difference operator Op
    % Op: (Difference Operator e.g. as received from GetStandardOp)
    % IdxCube: (tensor of order of problem dimension,
    %           size corresponding to problem discretization)
    %   e.g. for a simulation region (0,1)^3 with h1 = 0.1, h2=h3=0.2
    %   IdxCube must be a 10*5*5-tensor
    % rhs: (vector of length of numel(IdxCube))
    %   Right hand side of the PDE's finite difference linear equation
    %   system befor applying the BC
    % return value:
    % Op: modified difference operator
    % rhs: (vector of length of numel(IdxCube))
    %   modified right hand side of the PDE's finite difference linear
    %   equation system
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14
    Dim = length(size(IdxCube));

    if ~exist('rhs','var')
	rhs = zeros(numel(IdxCube),1);
    end

    str = char(' ' * ones(1,4*Dim-1));
    str(4*(1:Dim-1)) = ',';
    str(4*(1:Dim)-3) = ':';
    for iDim=1:Dim
        str((4*(iDim-1))+(1:3)) = ' 1 ';
        Idx = eval(['IdxCube(' str ')']);
        Op = ApplyStencil(IdxCube,1,Idx,Op);
        rhs(Idx) = 0;
        str((4*(iDim-1))+(1:3)) = 'end';
        Idx = eval(['IdxCube(' str ')']);
        Op = ApplyStencil(IdxCube,1,Idx,Op);
        rhs(Idx) = 0;
        str((4*(iDim-1))+(1:3)) = ' : ';
    end
end
