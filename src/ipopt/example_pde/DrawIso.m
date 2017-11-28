function DrawIso(varargin)
    % Draw surfaces of equal values of a function R^3->R
    % Possible interfaces:
    % DrawIso(v,x1,x2,x3,IsoVal)
    % DrawIso(v)
    % DrawIso(v,IsoVal)
    % DrawIso(v,x1,x2,x3)
    % DrawIso(v,x1,x2,x3,IsoVal) interpolates the 3D-matrix v at isovalues
    % Isoval (can be a vector) and plots it, where x1,x2,x3 are the
    % coordinates of the points in v. If no coordinates are given Indicees
    % are used instead. If no IsoVal is passed, the are set to
    % IsoVal = min(v(:)) + [1/3 2/3] * (max(v(:)) - min(v(:)));

    
    if ~ishold      % clear current axes if hold is off
        cla;
    end
    
    switch nargin
        case 1
            v = varargin{1};
        case 2
            v = varargin{1};
            IsoVal = varargin{2};
        case 4
            v = varargin{1};
            x1 = varargin{2};
            x2 = varargin{3};
            x3 = varargin{4};
        case 5
            v = varargin{1};
            x1 = varargin{2};
            x2 = varargin{3};
            x3 = varargin{4};
            IsoVal = varargin{5};
        otherwise
            error('wrong number of arguments')
    end

    if ~exist('x1','var')
        x1 = 1:size(v,1);
        x2 = 1:size(v,2);
        x3 = 1:size(v,3);
    end
    
    MinVal = min(v(:));
    MaxVal = max(v(:));
%    caxis([MinVal MaxVal]);
    ColMap = get(gcf,'Colormap');
    if ~exist('IsoVal','var')
        IsoVal = MinVal + [1/3 2/3]*(MaxVal-MinVal);
    end
    for iVal = 1:length(IsoVal)
        p = patch(isosurface(x2,x1,x3,v,IsoVal(iVal)));
        alpha(p,0.4);
        isonormals(x2,x1,x3,v,p);
        ColIdx = fix((IsoVal(iVal)-MinVal)/(MaxVal-MinVal)*size(ColMap,1))+1;
        ColIdx = min(size(ColMap,1),max(1,ColIdx));
        Color = ColMap(ColIdx,:);
        set(p,'FaceColor',Color,'EdgeColor','none');
    end
    daspect([1 1 1])
    camlight 
    % lighting gouraud
end