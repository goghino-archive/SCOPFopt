function rv = ApplyStencil(IdxCube,DiffStencil,PtIdx,Init)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Applies differnce stencil to points
    % IdxCube: (tensor of order of problem dimension,
    %           size corresponding to problem discretization)
    %   e.g. for a simulation region (0,1)^3 with h1 = 0.1, h2=h3=0.2
    %   IdxCube must be a 10*5*5-tensor
    % DiffStencil: (tensor of order of problem dimension,
    %           size corresponding to stencil size with reference point
    %           in the center)
    %   e.g. for a central difference approximation of d/dx1 in R^2 with h1=h2=0.1
    %   it looks DiffStencil = [0 -1 0; 0 0 0;0 1 0];
    % PtIdx: (vector) point indices (1-based) of points, where the stencil will be
    %   applied
    % Init: (Difference Operator e.g. as received from GetStandardOp or none)
    %   Difference operator that will be modified
    % return value:
    % rv: modified difference operator
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14    
    if ~exist('PtIdx','var')
		PtIdx = ':';
    end
	
    if ~exist('Init','var')
	Init = struct('DiagVal',{[]},'DiagIdx',{[]});
	rv = Init;
    else
	rv = Init;
	rv.DiagVal(PtIdx,:) = 0;
    end

    IdxNZDS = find(DiffStencil~=0);
    szDS = size(DiffStencil);
    if length(szDS)<ndims(IdxCube)
        szDS = [szDS ones(1,ndims(IdxCube)-length(szDS))];
    end
    SubNZ = Myind2sub(szDS,IdxNZDS);
    NumStencilEntries = size(SubNZ,1);
    if NumStencilEntries==0
	return
    end
    SubCenter = ceil(szDS/2);
    dSub = SubNZ-repmat(SubCenter,NumStencilEntries,1);
    szIdxCube = size(IdxCube);
    StencilDiags = Mysub2ind(szIdxCube,dSub)';
    NewDiagIdx = StencilDiags(~ismember(StencilDiags,Init.DiagIdx));
    NumNewDiags = length(NewDiagIdx);
    if NumNewDiags>0
	rv.DiagVal = [rv.DiagVal zeros(numel(IdxCube),NumNewDiags)];
	rv.DiagIdx = [rv.DiagIdx NewDiagIdx];
    end
    
    for iStencilEntry = 1:NumStencilEntries
	LogIdx = ( rv.DiagIdx==StencilDiags(iStencilEntry) );
	rv.DiagVal(PtIdx,LogIdx) = DiffStencil(IdxNZDS(iStencilEntry));
    end
end


function rv = Myind2sub(siz ,ndx)
	if isempty(ndx)
		rv = [];
		return
	end
	n = length(siz);
	rv = zeros(length(ndx),length(siz));
	k = [1 cumprod(siz(1:end-1))];

	for i = n:-1:1,
		vi = rem(ndx-1, k(i)) + 1;         
		vj = (ndx - vi)/k(i) + 1; 
		rv(:,i) = vj;
		ndx = vi;     
	end
end

function rv = Mysub2ind(sz,sub)
	rv = sub.*repmat([1 cumprod(sz(1:end-1))],size(sub,1),1);
	rv = sum(rv,2);
end
