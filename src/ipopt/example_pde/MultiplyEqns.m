function [Op, Rhs] = MultiplyEqns(Op, Coeff, Rhs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multiplies equations with values contained in Coeff
    % Op: (Difference Operator e.g. as received from GetStandardOp)
    %   Difference operator that will be multiplied with the discretized
    %   funtcion Coeff
    % Coeff: (vector of size of discretization mesh or single value)
    %   discretized function to multiply the discretized differential
    %   operator with
    % Rhs: (vector of size of discretization mesh or none)
    %   Right hand side of discretized PDE
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14

    if numel(Coeff) == 1
		Coeff = Coeff*ones(size(Op.DiagVal,1),1);
	end

	Op.DiagVal = Op.DiagVal.*repmat(Coeff,1,size(Op.DiagVal,2));
    if exist('Rhs')
        Rhs = Rhs.*repmat(Coeff,1,size(Rhs,2));
    end
end

