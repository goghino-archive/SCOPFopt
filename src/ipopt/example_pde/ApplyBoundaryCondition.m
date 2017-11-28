function [Op, Rhs] = ApplyBoundaryCondition(PtIdx,DiagIdxOut,NeumCoef,DiriCoef,g,Op,Rhs,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Applies Robin boundary condition (BC)
	% NeumCoef*dy/dn + DiriCoef y = g
    % to the difference operator Op and the
    % corresponding system right hand side Rhs, uses cerntral difference
    % quotient to approximate dy/dn
    % PtIdx: (vector) point indices (1-based) of points, where the BC will be
    %   applied
    % DiagIdxOut: (number) index offset of points, that should be eliminated by
    %   use of the BC, usually the difference of indices of the
    %   outlying point and its corresponding reference point
    % NeumCoef: (vector of length PtIdx or single number)
    %   Coefficient of normal derivative in BC
    % DiriCoef: (vector of length PtIdx or single number)
    %   Coefficient of function value in BC
    % g: (vector of length PtIdx or single number)
    %   right hand side in boundary condition
    % Op: (Difference Operator e.g. as received from GetStandardOp)
    % Rhs: (vector of length of mesh size or single number)
    %   Right hand side of the PDE's finite difference linear equation
    %   system befor applying the BC
    % h: (single value)
    %   discretization step
    % return value:
    % Op: modified difference operator
    % rhs: modified linear equation right hand side
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14
	if numel(NeumCoef)==1
		NeumCoef = NeumCoef*ones(numel(PtIdx),1);
	end
	if numel(DiriCoef)==1
		DiriCoef = DiriCoef*ones(numel(PtIdx),1);
	end
	if numel(g)==1
		g = g*ones(numel(PtIdx),1);
	end
	if numel(Rhs)==1
		Rhs = Rhs*ones(size(Op.DiagVal,1),1);
	end
	CenterDiag = find(Op.DiagIdx==0);
	if isempty(CenterDiag)
	    Op.DiagVal = [Op.DiagVal zeros(size(Op.DiagVal,1),1)];
	    Op.DiagIdx = [Op.DiagIdx 0];
	    CenterDiag = length(Op.DiagIdx);
	end
	%CenterDiag = (size(Diags,2)+1)/2;
	% First all boundary points, where NeumCoef = 0 (Dirichlet-BC)
	iA0 = NeumCoef==0;
	Idx = PtIdx(iA0);
	Op.DiagVal(Idx,:) = 0;
	Op.DiagVal(Idx,CenterDiag) = DiriCoef(iA0);
	Rhs(Idx,:) = repmat(g(iA0),1,size(Rhs,2));
	% All other points
	iAn0 = NeumCoef~=0;
	Idx = PtIdx(iAn0);
	if ~isempty(Idx)
		OuterDiag = find(Op.DiagIdx==DiagIdxOut);
		if ~isempty(OuterDiag)
			InnerDiag = find(Op.DiagIdx==-DiagIdxOut);
			if isempty(InnerDiag)
			    Op.DiagVal = [Op.DiagVal zeros(size(Op.DiagVal,1),1)];
			    Op.DiagIdx = [Op.DiagIdx -DiagIdxOut];
			    InnerDiag = length(DiagIdx);
			end
			Op.DiagVal(Idx,InnerDiag) = Op.DiagVal(Idx,InnerDiag) + Op.DiagVal(Idx,OuterDiag);
			Op.DiagVal(Idx,CenterDiag) = Op.DiagVal(Idx,CenterDiag) - 2*h*Op.DiagVal(Idx,OuterDiag).*DiriCoef(iAn0)./NeumCoef(iAn0);
			Rhs(Idx,:) = Rhs(Idx,:) - 2*h*repmat(Op.DiagVal(Idx,OuterDiag).*g(iAn0)./NeumCoef(iAn0),1,size(Rhs,2));
			Op.DiagVal(Idx,OuterDiag) = 0;
		end
	end
end
