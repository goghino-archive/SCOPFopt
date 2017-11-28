function A = CreateSystemMatrix(Op)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creates a system matrix from the finite difference Operator, such,
    % that applying this finite difference operator if matrix
    % multiplication
    % Op: (Difference Operator e.g. as received from GetStandardOp)
    %   Difference operator that will be used to create the system matrix
    % return value:
    % A: system matrix
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14    
    A = spdiags(Op.DiagVal,-Op.DiagIdx,size(Op.DiagVal,1),size(Op.DiagVal,1)).';
end

