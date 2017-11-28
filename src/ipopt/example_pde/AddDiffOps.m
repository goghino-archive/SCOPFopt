function DSDRv = AddDiffOps(DSD1,DSD2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Adds difference operators DSD1 and DSD2
    % DSD1: Difference Operator1 e.g. from GetStandardOp
    % DSD2: Difference Operator2 e.g. from GetStandardOp
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14
    DSDRv = struct('DiagVal',{[]},'DiagIdx',{[]});
    DSDRv.DiagIdx = union(DSD1.DiagIdx,DSD2.DiagIdx);
	DSDRv.DiagVal = zeros(size(DSD1.DiagVal,1),numel(DSDRv.DiagIdx));
	for iDiag = 1:numel(DSDRv.DiagIdx)
		iD1 = find(DSD1.DiagIdx==DSDRv.DiagIdx(iDiag));
		iD2 = find(DSD2.DiagIdx==DSDRv.DiagIdx(iDiag));
		% it is not possible, that iD1 and iD2 is empty, since DiagIdx = union(DiagIdx1,DiagsIdx2);
		if isempty(iD1)
			DSDRv.DiagVal(:,iDiag) = DSD2.DiagVal(:,iD2);
		elseif isempty(iD2)
			DSDRv.DiagVal(:,iDiag) = DSD1.DiagVal(:,iD1);
		else
			DSDRv.DiagVal(:,iDiag) = DSD1.DiagVal(:,iD1)+DSD2.DiagVal(:,iD2);
		end
	end
end
