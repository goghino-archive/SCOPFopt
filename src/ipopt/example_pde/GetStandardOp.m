function rv = GetStandardOp(strID,IdxCube,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % returns a finite difference approximation for several well known
    % differential operators
    % see also 'GetStandardStencil'
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14    
	Stencil = GetStandardStencil(strID,IdxCube,h);
	rv = ApplyStencil(IdxCube,Stencil);
end
