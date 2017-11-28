function rv = GetStandardStencil(strID,IdxCube,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % returns a finite difference approximation stencil 
    % for several well known differential operators
    % strID: string
    %   'CentralDxi': Central difference quotient for approximation of
    %       d/dxi
    %   'OneSideDxi+': Right side difference quotient for approximation of
    %       d/dxi
    %   'OneSideDxi-': Left side difference quotient for approximation of
    %       d/dxi
    %   'OneSideDxi++': Right side difference quotient for approximation of
    %       d/dxi, with constistency order 2
    %   'OneSideDxi--': Left side difference quotient for approximation of
    %       d/dxi, with constistency order 2
    %   'Laplace': Standard difference stencil for Laplacian of constistency order 2
    % All stencils are created with respect to the problem dimension
    % according to IdxCube
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_14    	
    Dim = length(size(IdxCube));
	if length(h)==1 && Dim>1
		h = repmat(h,1,Dim);
    end
    strID = upper(strID);
    
    [Val, count, errmsg] = sscanf(strID,'CENTRALDX%d');
    if count==1 && strcmp(errmsg,'')
        iDxDim = Val(1);
	if iDxDim>Dim
	    error('dimension exceeds index cube');
	end
        strID = 'CENTRALDX';
    end

    [Val, count, errmsg] = sscanf(strID,'ONESIDEDX%d%[+-]');
    if count>=2 && strcmp(errmsg,'')
        iDxDim = Val(1);
	if iDxDim>Dim
	    error('dimension exceeds index cube');
	end
        Dir = char(Val(2:end)');
        strID = 'ONESIDEDX';
    end
    
    [Val, count, errmsg] = sscanf(strID,'D2X%d');
    if count==1 && strcmp(errmsg,'')
        iDxDim = Val(1);
	if iDxDim>Dim
	    error('dimension exceeds index cube');
	end
        strID = 'D2X';
    end

    Stencil = [];
	switch(strID)
		case 'D2X'
			str = char(' ' * ones(1,2*Dim-1));
			str(2*(1:Dim-1)) = ',';
			str(2*(1:Dim)-1) = '3';
			Stencil = eval(['zeros(' str ')']);
			str(2*(1:Dim)-1) = '2';
			str(2*iDxDim-1) = '1';
			eval(['Stencil(' str ') = 1/(h(iDxDim)^2);']);
			str(2*iDxDim-1) = '3';
			eval(['Stencil(' str ') = 1/(h(iDxDim)^2);']);
			str(2*iDxDim-1) = '2';
			eval(['Stencil(' str ') = -2/(h(iDxDim)^2);']);
		case 'LAPLACE'
			str = char(' ' * ones(1,2*Dim-1));
			str(2*(1:Dim-1)) = ',';
			str(2*(1:Dim)-1) = '3';
			Stencil = eval(['zeros(' str ')']);
			str(2*(1:Dim)-1) = '2';
			for iDim=1:Dim
				str(2*iDim-1) = '1';
				eval(['Stencil(' str ') = 1/h(iDim)^2;']);
				str(2*iDim-1) = '3';
				eval(['Stencil(' str ') = 1/h(iDim)^2;']);
				str(2*iDim-1) = '2';
				eval(['Stencil(' str ') = Stencil(' str ') - 2/h(iDim)^2;']);
			end
		case 'CENTRALDX'
			str = char(' ' * ones(1,2*Dim-1));
			str(2*(1:Dim-1)) = ',';
			str(2*(1:Dim)-1) = '3';
			Stencil = eval(['zeros(' str ')']);
			str(2*(1:Dim)-1) = '2';
			str(2*iDxDim-1) = '1';
			eval(['Stencil(' str ') = -1/(2*h(iDxDim));']);
			str(2*iDxDim-1) = '3';
			eval(['Stencil(' str ') = 1/(2*h(iDxDim));']);
		case 'ONESIDEDX'
			switch(Dir)
				case '+' % simple foreward difference quotient
				    str = char(' ' * ones(1,2*Dim-1));
				    str(2*(1:Dim-1)) = ',';
				    str(2*(1:Dim)-1) = '3';
				    Stencil = eval(['zeros(' str ')']);
				    str(2*(1:Dim)-1) = '2';
				    eval(['Stencil(' str ') = -1/h(iDxDim);']);
				    str(2*iDxDim-1) = '3';
				    eval(['Stencil(' str ') =  1/h(iDxDim);']);
				case '-'  % simple backward difference quotient
				    str = char(' ' * ones(1,2*Dim-1));
				    str(2*(1:Dim-1)) = ',';
				    str(2*(1:Dim)-1) = '3';
				    Stencil = eval(['zeros(' str ')']);
				    str(2*(1:Dim)-1) = '2';
				    eval(['Stencil(' str ') =  1/h(iDxDim);']);
				    str(2*iDxDim-1) = '1';
				    eval(['Stencil(' str ') = -1/h(iDxDim);']);
				case '++'
				    str = char(' ' * ones(1,2*Dim-1));
				    str(2*(1:Dim-1)) = ',';
				    str(2*(1:Dim)-1) = '5';
				    Stencil = eval(['zeros(' str ')']);
				    str(2*(1:Dim)-1) = '3';
				    eval(['Stencil(' str ') = -3/(2*h(iDxDim));']);
				    str(2*iDxDim-1) = '4';
				    eval(['Stencil(' str ') =  4/(2*h(iDxDim));']);
				    str(2*iDxDim-1) = '5';
				    eval(['Stencil(' str ') = -1/(2*h(iDxDim));']);
				case '--'
				    str = char(' ' * ones(1,2*Dim-1));
				    str(2*(1:Dim-1)) = ',';
				    str(2*(1:Dim)-1) = '5';
				    Stencil = eval(['zeros(' str ')']);
				    str(2*(1:Dim)-1) = '3';
				    eval(['Stencil(' str ') = 3/(2*h(iDxDim));']);
				    str(2*iDxDim-1) = '2';
				    eval(['Stencil(' str ') = -4/(2*h(iDxDim));']);
				    str(2*iDxDim-1) = '1';
				    eval(['Stencil(' str ') = 1/(2*h(iDxDim));']);
				otherwise
				    error('unknown difference operator')
			end
		case '1'
			str = char(' ' * ones(1,2*Dim-1));
			str(2*(1:Dim-1)) = ',';
			str(2*(1:Dim)-1) = '3';
			Stencil = eval(['zeros(' str ')']);
			str(2*(1:Dim)-1) = '2';
			eval(['Stencil(' str ') = 1;']);
		case '0'
			str = char(' ' * ones(1,2*Dim-1));
			str(2*(1:Dim-1)) = ',';
			str(2*(1:Dim)-1) = '1';
			Stencil = eval(['zeros(' str ')']);
		otherwise
		    error('unknown difference operator');
    end
    rv = Stencil;
end
