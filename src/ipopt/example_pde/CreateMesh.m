function varargout = CreateMesh(Data)
    sz = size(Data);
    if sz(2) ~=3
	error('Data must have 3 columns');
    end
    if nargout~=sz(1)+1
        error(['For a ' num2str(sz(1)) 'D mesh there must be ' num2str(sz(1)+1) ' return values']);
    end
    str = ['['];
    for iDim = 1:sz(1)
        str = [str 'x' num2str(iDim) ','];
    end
    str(end) = ']';
    
    str = [str '=ndgrid('];
    for iDim = 1:sz(1)
        str = [str 'Data(' num2str(iDim) ',1):Data(' num2str(iDim) ',2):Data(' num2str(iDim) ',3),'];
    end
    str(end) = ')';
    str = [str ';'];
    eval(str);

    varargout = cell(1,nargout);
    for iDim = 1:sz(1)
        varargout(iDim) = {eval(['x' num2str(iDim)])};
    end
    
    IdxCube = reshape(1:numel(x1),size(x1));
    varargout{sz(1)+1} = IdxCube;
end
