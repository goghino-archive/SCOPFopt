function SliceGUI(Val,x1,x2,x3,SlPosx1,SlPosx2,SlPosx3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interactive visualization tool of a function R^3->R
    % ATTENTION: values are asumed to be discretized from a mesh created by
    % ndgrid, NOT by meshgrid, i.e. different mesh point ordering
    % Val: (tensor of order 3)
    %   function values
    % x1,x2,x3: (tensor of order 3 or none)
    % grid coordinates as received by ndgrid
    %   Author: Johannes Huber
    %           Mathematische Institut
    %           University of Basel
    %   2009_01_29
    if ~exist('x1','var')
        D1 = 1:size(Val,1);
        D2 = 1:size(Val,2);
        D3 = 1:size(Val,3);
        [x1,x2,x3] = ndgrid(D1,D2,D3);
    end
    Sz = size(Val);

    % Change Dimension order to make change ndgrid -> meshgrid
    x1 = permute(x1,[2 1 3]);
    x2 = permute(x2,[2 1 3]);
    x3 = permute(x3,[2 1 3]);
    Val = permute(Val,[2 1 3]);
    
    Min = min(Val(:));
    Max = max(Val(:));
    if Max<=Min
        Max = Min+1;
    end

    if exist('SlPosx1','var')
        SlPosx1 = min(SlPosx1,max(x1(1,:,1)));
        SlPosx1 = max(SlPosx1,min(x1(1,:,1)));
        Idx = find((abs(x1(1,:,1)-SlPosx1)-eps)<=(x1(1,2,1)-x1(1,1,1))/2);
        SlPosx1 = Idx(1);
    else
        SlPosx1 = ceil(Sz(1)/2);
    end
    
    if exist('SlPosx2','var')
        SlPosx2 = min(SlPosx2,max(x2(:,1,1)));
        SlPosx2 = max(SlPosx2,min(x2(:,1,1)));
        Idx = find((abs(x2(:,1,1)-SlPosx2)-eps)<=(x2(2,1,1)-x2(1,1,1))/2);
        SlPosx2 = Idx(1);
    else
        SlPosx2 = ceil(Sz(2)/2);
    end

    if exist('SlPosx3','var')
        SlPosx3 = min(SlPosx3,max(x3(1,1,:)));
        SlPosx3 = max(SlPosx3,min(x3(1,1,:)));
        Idx = find((abs(x3(1,1,:)-SlPosx3)-eps)<=(x3(1,1,2)-x3(1,1,1))/2);
        SlPosx3 = Idx(1);
    else
        SlPosx3 = ceil(Sz(3)/2);
    end
    
    figure('Toolbar','figure');
    hMainAx = axes('OuterPosition',[0    0   0.6  1]);

    yOffset = 0.05;
    hX3PlaneAx = axes('OuterPosition',[0.75 yOffset 0.2 0.25]);
    hTextX3 = annotation('textbox','edgecolor','none',...
        'units','normal','Position',[0.65,yOffset+0.22,0.1,0.05],...
        'String',['x3: ' num2str(x3(1,1,SlPosx3))]);
    hSliderX3 = uicontrol('Style','slider',...
        'units','normal','Position',[0.65,yOffset,0.05,0.2],...
        'Min',0,'Max',Sz(3)+1,...
        'SliderStep',[1/(Sz(3)+1) 1/(Sz(3)+1)],'Value',SlPosx3,...
        'Callback',@OnSliderChanged);

    yOffset = 0.37;
    hX2PlaneAx = axes('OuterPosition',[0.75 yOffset 0.2 0.25]);
    hTextX2 = annotation('textbox','edgecolor','none',...
        'units','normal','Position',[0.65,yOffset+0.22,0.1,0.05],...
        'String',['x2: ' num2str(x2(SlPosx2,1,1))]);
    hSliderX2 = uicontrol('Style','slider',...
        'units','normal','Position',[0.65,yOffset,0.05,0.2],...
        'Min',0,'Max',Sz(2)+1,'SliderStep',[1/(Sz(2)+1) 1/(Sz(2)+1)],'Value',SlPosx2,...
        'Callback',@OnSliderChanged);

    yOffset = 0.68;
    hX1PlaneAx = axes('OuterPosition',[0.75 yOffset 0.2 0.25]);
    hTextX1 = annotation('textbox','edgecolor','none',...
        'units','normal','Position',[0.65,yOffset+0.22,0.1,0.05],...
        'String',['x1: ' num2str(x1(1,SlPosx1,1))]);
    hSliderX1 = uicontrol('Style','slider',...
        'units','normal','Position',[0.65,yOffset,0.05,0.2],...
        'Min',0,'Max',Sz(1)+1,'SliderStep',[1/(Sz(1)+1) 1/(Sz(1)+1)],'Value',SlPosx1,...
        'Callback',@OnSliderChanged);

     DrawSlices(1);

    function OnSliderChanged(hObject, eventdata)
        SlPosx1 = round(get(hSliderX1,'Value'));
        if SlPosx1<1 || SlPosx1>Sz(1)
            SlPosx1 = [];
        end
        set(hTextX1,'String',['x1: ' num2str(x1(1,SlPosx1,1))]);

        SlPosx2 = round(get(hSliderX2,'Value'));
        if SlPosx2<1 || SlPosx2>Sz(2)
            SlPosx2 = [];
        end
        set(hTextX2,'String',['x2: ' num2str(x2(SlPosx2,1,1))]);

        SlPosx3 = round(get(hSliderX3,'Value'));
        if SlPosx3<1 || SlPosx3>Sz(3)
            SlPosx3 = [];
        end
        set(hTextX3,'String',['X3: ' num2str(x3(1,1,SlPosx3))]);
        DrawSlices(0);
    end
    
    function DrawSlices(Init)
        axes(hMainAx);
        if(~Init)
            CP = campos;
        end
        slice(x1,x2,x3,Val,x1(1,SlPosx1,1),x2(SlPosx2,1,1),x3(1,1,SlPosx3));    % Draw the slices
        if(~Init)
            campos(CP);
        end
        axis([x1(1,1,1) x1(end,end,end) x2(1,1,1) x2(end,end,end) x3(1,1,1) x3(end,end,end) Min Max]);
        xlabel('x1');
        ylabel('x2');
        zlabel('x3');
        colorbar;
        
        % find slice planes' handles
        tmp = get(hMainAx,'Children');
        X1Plane = [];
        X2Plane = [];
        X3Plane = [];
        for ih = tmp'
            XData = get(ih,'XData');
            if XData(1,1)==XData(1,2) && XData(1,1)==XData(2,1)
                X1Plane = ih;
            else
                YData = get(ih,'YData');
                if YData(1,1)==YData(1,2) && YData(1,1)==YData(2,1)
                    X2Plane = ih;
                else
                    ZData = get(ih,'ZData');
                    if ZData(1,1)==ZData(1,2) && ZData(1,1)==ZData(2,1)
                        X3Plane = ih;
                    else
                        warning('Unknown Child');
                    end
                end
            end
        end

        % plot small figures
        if ~isempty(X1Plane)
            if ishandle(X1Plane)
                tmp = get(X1Plane,'CData');
                axes(hX1PlaneAx);
                x2Plot = squeeze(x2(:,1,:));
                x3Plot = squeeze(x3(:,1,:));
                mesh(x2Plot,x3Plot,tmp);
                zlim([Min Max]);
                xlabel('x2');
                ylabel('x3');
            end
        else
            axes(hX1PlaneAx);
            cla;
        end
        
        if ~isempty(X2Plane)
            if ishandle(X2Plane)
                tmp = get(X2Plane,'CData');
                axes(hX2PlaneAx);
                x1Plot = squeeze(x1(1,:,:));
                x3Plot = squeeze(x3(1,:,:));
                mesh(x3Plot,x1Plot,tmp);
                xlabel('x3');
                ylabel('x1');
                zlim([Min Max]);
            end
        else
            axes(hX2PlaneAx);
            cla;
        end
        if ~isempty(X3Plane)
            if ishandle(X3Plane)
                tmp = get(X3Plane,'CData');
                axes(hX3PlaneAx);
                x1Plot = squeeze(x1(:,:,1));
                x2Plot = squeeze(x2(:,:,1));
                mesh(x1Plot,x2Plot,tmp);
                xlabel('x1');
                ylabel('x2');
                zlim([Min Max]);
            end
        else
            axes(hX3PlaneAx);
            cla;
        end
        axes(hMainAx);  % setting axes properties should go to the main axes, so activate it after drawing (e.g. title(.) in workspace
        drawnow;
    end
end
