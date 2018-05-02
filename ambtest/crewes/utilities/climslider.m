function climslider(haxe,hfig,pos,N,xn)
%
% haxe ... handle of the axes to control
% hfig ... handle of the figure to put the tool in
% pos ... position in hfig (normalized) of the tool
% N,xn ... the return values from [N,xn]=hist(data(:),100)
%
% Refresh the histogram
% climslider('refresh',hfig,N,xn)
% here hfig is the handle of the Figure containing the slider, N and xn are the new histogram

global DRAGLINE_MOTION DRAGLINE_XLIMS DRAGLINE_YLIMS DRAGLINE_SHOWPOSN DRAGLINE_CALLBACK DRAGLINE_MOTIONCALLBACK DRAGLINE_PAIRED 

if(ischar(haxe))
    action=haxe;
else
    action='init';
end


if(strcmp(action,'init'))
    hclim=findobj(hfig,'tag','clim');
    if(isempty(hclim))
        if(verLessThan('matlab','R2017a'))
            figure(hfig);
            hax=axes('position',[.1 .2 .8 .6]);
            hclim=uicontrol(hfig,'style','text','visible','off','tag','clim');
        else
            hclim=uipanel(hfig,'position',pos,'tag','clim');
            hax=axes(hclim,'position',[.1 .2 .8 .6]);
        end
        kol2=[1 0 0];%color of clim lines
        clim=get(haxe,'clim');
        %horizontal orientation
        bar(xn,N);
        set(hax,'ytick',[]);ylabel('histogram')
        xl=get(gca,'xlim');
        xrange=xl(2)-xl(1);
        yl=get(gca,'ylim');
        line([clim(1) clim(1)],yl,'color',kol2,'buttondownfcn','climslider(''dragline'');',...
            'tag','clim1');
        line([clim(2) clim(2)],yl,'color',kol2,'buttondownfcn','climslider(''dragline'');',...
            'tag','clim2');
        set(hax,'xlim',[min([xl(1)-.2*xrange, clim(1)-.2*xrange]) max([xl(2)+.2*xrange clim(2)+.2*xrange])]);
        title('drag the red lines');
        set(hclim,'userdata',haxe);
    end
    hax=findobj(hfig,'type','axes');
    axes(hax);
    climslider('setclim');
        
elseif(strcmp(action,'setclim'))
    h1=findobj(gca,'tag','clim1');
    %hclim=get(gca,'parent');
    hclim=findobj(gcf,'tag','clim');
    haxe=get(hclim,'userdata');
    xx=get(h1,'xdata');
    clim1=xx(1);
%     yy=get(h1,'ydata');
%     if(diff(yy)==1)
%         clim1=xx(1);
%     else
%         clim1=yy(1);
%     end
    h2=findobj(gca,'tag','clim2');
    xx=get(h2,'xdata');
    clim2=xx(2);
%     yy=get(h2,'ydata');
%     if(diff(yy)==1)
%         clim2=xx(1);
%     else
%         clim2=yy(1);
%     end
    set(haxe,'clim',[clim1 clim2])
elseif(strcmp(action,'dragline'))
    hh=gco;
    xl=get(gca,'xlim');
    h1=findobj(gca,'tag','clim1');
    h2=findobj(gca,'tag','clim2');
    xx=get(h1,'xdata');
    clim1=xx(1);
%     yy=get(h1,'ydata');
%     if(diff(yy)==1)
%         clim1=xx(1);
%     else
%         clim1=yy(1);
%     end
    xx=get(h2,'xdata');
    clim2=xx(2);
%     yy=get(h2,'ydata');
%     if(diff(yy)==1)
%         clim2=xx(1);
%     else
%         clim2=yy(1);
%     end
    DRAGLINE_MOTION='xonly';
    if(hh==h1)
        %we are dragging h1
        DRAGLINE_XLIMS=[xl(1) clim2];
        DRAGLINE_PAIRED=h2;
    else
        %we are dragging h2
        DRAGLINE_XLIMS=[clim1 xl(2)];
        DRAGLINE_PAIRED=h1;
    end
    DRAGLINE_YLIMS=[];
    DRAGLINE_SHOWPOSN='on';
    DRAGLINE_CALLBACK='';
    DRAGLINE_MOTIONCALLBACK='climslider(''setclim'');';
    dragline('click')
elseif(strcmp(action,'refresh'))
    xn=N;
    N=pos;
    figure(hfig);
    haxe=findobj(hfig,'type','axes');
    hbar=findobj(haxe,'type','bar');
    delete(hbar);
    set(hfig,'currentaxes',haxe);
    hold on
    bar(xn,N);
    hold off
end