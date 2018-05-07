function datar=seisplotfdom(seis,t,x,twin,tinc,fmt0,dname,cmap)
% seisplotfdom: examine time-variant dominant frequency
%
% datar=seisplotfdom(seis,t,x,twin,tinc,fmt0,dname,cmap)
%
% A new figure is created and divided into two axes (side-by-side). The first axes shows the
% input seismic gather and the second shows the dominant frequency section.
%
% seis ... 2D seismic matrix that spectral decomp was done on
% seissd ... 3D seismic matrix that resulted from spectral decomp. Can be either amplitude or phase.
% t ... time coordinate vector for seis
% x ... space coordinate vector for seis
% twin ... width (seconds) of the Gaussian window (standard deviation)
% tinc ... temporal shift (seconds) between windows
% fmt0 ... length 2 vector specifying a maximum signal frequency at some
%           time. For example, if signal is known to be bounded by 100Hz at 1 second then
%           fmt0=[100,1]. This determines the maximum frequency of integration at that time. For all
%           other times, fm is hyperbolically interpolated. That is, for time tk, fmk will be
%           fmk=fmt0(1)*fmt0(2)/tk.
%           NOTE: fmt0 can be specified as a single scalar in which case it is assumed the maximum
%           frequency is fmt0 and not time variant.
% dname ... text string nameing the first seismic matrix. Enter [] or '' for no name.
% cmap ... name of starting colormap to use
% ********** default 'seisclrs' *********
%
% datar ... Return data which is a length 2 cell array containing
%           datar{1} ... handle of the first seismic axes
%           datar{2} ... handle of the fdom axis
% These return data are provided to simplify plotting additional lines and
% text.
% 
% G.F. Margrave, Devon, 2018
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

%global DRAGLINE_MOTION DRAGLINE_XLIMS DRAGLINE_YLIMS DRAGLINE_SHOWPOSN DRAGLINE_CALLBACK DRAGLINE_MOTIONCALLBACK DRAGLINE_PAIRED %#ok<NUSED>
global NEWFIGVIS
if(~ischar(seis))
    action='init';
else
    action=seis;
end

datar=[];%initialize return data to null

if(strcmp(action,'init'))
    
    if(length(t)~=size(seis,1))
        error('time coordinate vector does not match first seismic matrix');
    end
    if(length(x)~=size(seis,2))
        error('space coordinate vector does not match first seismic matrix');
    end
    if(length(fmt0)~=2)
        error('fmt0 must be of length 2');
    end
    if(~between(t(1),t(end),fmt0(2),2))
        error('fmt0(2) must lie between t(1) and t(end)')
    end
    dt=t(2)-t(1);
    fnyq=.5/dt;
    fmax=.5*fnyq;
    tfmax=mean(t);
    if(~between(0,fnyq,fmt0(1)))
        error('fmt0(1) must lie between 0 and Nyquist')
    end
    
    if(nargin<7)
        dname=[];
    end
    if(nargin<8)
        cmap='seisclrs';
    end

    xwid=.35;
    xwid2=.35;
    yht=.8;
    xsep=.05;
    xnot=.5*(1-xwid-xwid2-1.5*xsep);
    ynot=.1;
    
    if(~isempty(NEWFIGVIS))
        figure('visible',NEWFIGVIS);
    else
        figure
    end
    hax1=subplot('position',[xnot ynot xwid yht]);

    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seis);
    clim=[am-clip*sigma am+clip*sigma];
        
    imagesc(x,t,seis,clim);
    brighten(.5);
    grid
    title(dname,'interpreter','none')
    maxmeters=7000;
    if(max(t)<10)
        ylabel('time (s)')
    elseif(max(t)<maxmeters)
        ylabel('depth (m)')
    else
        ylabel('depth (ft)')
    end
    if(max(x)<maxmeters)
        xlabel('distance (m)')
    else
        xlabel('distance (ft)')
    end
    %make a clip control

    xnow=xnot+xwid;
    wid=.055;ht=.05;
    ynow=ynot+yht;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clip1','units','normalized',...
        'position',[xnow,ynow,.8*wid,ht],'callback','seisplotfdom(''clip1'');','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,hax1},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    %make a help button
    uicontrol(gcf,'style','pushbutton','string','Info','tag','info','units','normalized',...
        'position',[xnot,ynow+ht,.5*wid,.5*ht],'callback','seisplotfdom(''info'');',...
        'backgroundcolor','y');
    
    %the hide seismic button
    xnow=xnot;
    uicontrol(gcf,'style','pushbutton','string','Hide seismic','tag','hideshow','units','normalized',...
        'position',[xnow,ynow,wid,.5*ht],'callback','seisplotfdom(''hideshow'');','userdata','hide');
    %the toggle button
    ynow=ynow+.5*ht;
    uicontrol(gcf,'style','pushbutton','string','Toggle','tag','toggle','units','normalized',...
        'position',[xnow,ynow,wid,.5*ht],'callback','seisplotfdom(''toggle'');','visible','off');
    
%     ht=.5*ht;
%     ynow=ynow-sep;
%     uicontrol(gcf,'style','pushbutton','string','brighten','tag','brighten','units','normalized',...
%         'position',[xnow,ynow,wid,ht],'callback','seisplotfdom(''brighten'')',...
%         'tooltipstring','push once or multiple times to brighten the images');
%     ynow=ynow-ht-sep;
%     uicontrol(gcf,'style','pushbutton','string','darken','tag','darken','units','normalized',...
%         'position',[xnow,ynow,wid,ht],'callback','seisplotfdom(''brighten'')',...
%         'tooltipstring','push once or multiple times to darken the images');
%     ynow=ynow-ht-sep;
%     uicontrol(gcf,'style','text','string','lvl 0','tag','brightness','units','normalized',...
%         'position',[xnow,ynow,wid,ht],...
%         'tooltipstring','image brightness (both images)','userdata',0);
    
    set(hax1,'tag','seis');
    
    hax2=subplot('position',[xnot+xwid+xsep ynot xwid2 yht]);
    set(hax2,'tag','seisfd');
%     [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seissd);
    
    
    
    %make a clip control
    sep=.01;
    xnow=xnot+2*xwid+xsep+sep+wid;
    ht=.05;
    ynow=ynot+yht;
    %wid=.045;sep=.005;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clip2','units','normalized',...
        'position',[xnow,ynow,.8*wid,ht],'callback','seisplotfdom(''clip2'');','value',iclip,...
        'userdata',{clips,zeros(1,4)},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    ht=.025;
    ynow=ynow-.5*ht;
    %colormap controls
    uicontrol(gcf,'style','text','string','Colormap:','tag','colomaplabel',...
        'units','normalized','position',[xnow ynow wid ht]);
    if(exist('parula','file')==2)
        colormaps={'seisclrs','redblue','redblue2','redblue3','blueblack','bluebrown','greenblack',...
            'greenblue','jet','parula','copper','bone','gray','winter'};
    else
        colormaps={'seisclrs','redblue','redblue2','redblue3','blueblack','bluebrown','greenblack',...
            'greenblue','jet','copper','bone','gray','winter'};
    end
    icolor=0;
    for k=1:length(colormaps)
        if(strcmp(colormaps{k},cmap))
            icolor=k;
        end
    end
    if(icolor==0)
        nk=length(colormaps);
        colormaps{nk+1}=cmap;
        icolor=nk+1;
    end
    ynow=ynow-ht;
    uicontrol(gcf,'style','popupmenu','string',colormaps,'tag','colormap',...
        'units','normalized','position',[xnow ynow wid ht],'callback',...
        'seisplotfdom(''colormap'');','value',icolor);
    ynow=ynow-2*ht-2.5*sep;
    hbg=uibuttongroup('position',[xnow,ynow,wid,3*ht],'title','Colormap goes to','tag','cmapgt');
    uicontrol(hbg,'style','radiobutton','string','left','tag','left','units','normalized',...
        'position',[0 2/3 1 1/3],'value',0);
    uicontrol(hbg,'style','radiobutton','string','right','tag','right','units','normalized',...
        'position',[0 1/3 1 1/3],'value',1);
    uicontrol(hbg,'style','radiobutton','string','both','tag','both','units','normalized',...
        'position',[0 0 1 1/3],'value',0);
    
    %controls to choose the dominant frequency display section 
    ynow=ynow-2*ht-4*sep;
    hbg=uibuttongroup('position',[xnow,ynow,wid,3*ht],'title','Display choice','tag','choices',...
        'selectionchangedfcn','seisplotfdom(''choice'');');
    ww=1;
    hh=.333;
    uicontrol(hbg,'style','radiobutton','string','Frequency','units','normalized','tag','freq',...
        'position',[0,2*hh,ww,hh],'value',1,'tooltipstring','Display dominant frequency');
    uicontrol(hbg,'style','radiobutton','string','Bandwidth','units','normalized','tag','bw',...
        'position',[0,hh,ww,hh],'value',0,'tooltipstring','Display bandwidth about dominant frequency');
    uicontrol(hbg,'style','radiobutton','string','Amplitude','units','normalized','tag','amp',...
        'position',[0,0,ww,hh],'value',0,'tooltipstring','Display amplitude at dominant frequency');
    
    %dominant frequency parameters
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','Fdom parameters:','units','normalized',...
        'position',[xnow,ynow,wid,ht]);
    ynow=ynow-ht-sep;
    wid=wid*.4;
    uicontrol(gcf,'style','text','string','Twin:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring',...
        'The standard deviation of the Gaussian window (seconds)');
    uicontrol(gcf,'style','edit','string',num2str(twin),'units','normalized','tag','twin',...
        'position',[xnow+wid+sep,ynow,wid,ht],'tooltipstring','Enter a value in seconds around 5 to 10 times the sample rate');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','Tinc:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring',...
        'The tewmporal increment between Gaussian windows (seconds)');
    uicontrol(gcf,'style','edit','string',num2str(tinc),'units','normalized','tag','tinc',...
        'position',[xnow+wid+sep,ynow,wid,ht],'tooltipstring','Enter a value smaller than Twin but not smaller than the sample rate');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','Fmax:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring',...
        'Should be a value in Hz about 50% larger than the maximum signal frequency');
    uicontrol(gcf,'style','edit','string',num2str(fmax),'units','normalized','tag','fmax',...
        'position',[xnow+wid+sep,ynow,wid,ht],'tooltipstring','Frequencies larger than this are ignored');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','Tfmax:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring',...
        'Time at which Fmax is specified');
    uicontrol(gcf,'style','edit','string',num2str(tfmax),'units','normalized','tag','tfmax',...
        'position',[xnow+wid+sep,ynow,wid,ht],'tooltipstring','If in doubt, don''t touch');
    ynow=ynow-ht-sep;
    wid=0.055;
    uicontrol(gcf,'style','pushbutton','string','Apply','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotfdom(''apply'');',...
        'tooltipstring','Apply current fdom specs');
    
    ynow=ynow-2*ht-sep;
     uicontrol(gcf,'style','text','string','Compute performace:','units','normalized',...
        'position',[xnow,ynow,1.2*wid,ht]);
    ynow=ynow-ht-sep;
     uicontrol(gcf,'style','text','string','','units','normalized','tag','performance',...
        'position',[xnow,ynow,1.2*wid,ht]);
    
    
    %results popup
    pos=get(hax2,'position');
    xnow=pos(1);
    ynow=pos(2)+pos(4)-ht;
    wid=pos(3);
    ht=3*ht;
    fs=16;
    uicontrol(gcf,'style','popupmenu','string','Diddley','units','normalized','tag','results',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotfdom(''select'');','fontsize',fs,...
        'fontweight','bold');
    
    %zoom buttons
%     wid=.1;
%     pos=get(hax,'position');
%     xnow=pos(1)+.5*pos(3)-.5*wid;
%     ynow=.97;
%     uicontrol(gcf,'style','pushbutton','string','Zoom #1 like #2','units','normalized',...
%         'position',[xnow ynow wid .5*ht],'tag','1like2','callback','seisplotfdom(''equalzoom'');');
%     
%     pos=get(hax,'position');
%     xnow=pos(1)+.5*pos(3)-.5*wid;
%     uicontrol(gcf,'style','pushbutton','string','Zoom #2 like #1','units','normalized',...
%         'position',[xnow ynow wid .5*ht],'tag','2like1','callback','seisplotfdom(''equalzoom'');');
    seisplotfdom('apply');
    colormap(seisclrs);
    brighten(.5)
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
    
   
    set(gcf,'name',['Dominant frequency for ' dname],'closerequestfcn','seisplotfdom(''close'');',...
        'userdata',-999.25);%userdata here is just a placeholder to cause plotimage3D to do the right thing
    if(nargout>0)
        datar=cell(1,2);
        datar{1}=hax1;
        datar{2}=hax2;
    end
elseif(strcmp(action,'info'))
    hthisfig=gcf;
    msg=['The axes at left (the seismic axes) shows the ordinary sesimic, the axes at right ',...
        '(fdom axes) shows one of the three dominant frequency sections. At far left above the seismic axes is a ',...
        'button labelled "Hide seismic". Clicking this removes the seismic axes from the display ',...
        'allows the fdom axes to fill the window. This action also displays a new button labelled ',...
        '"Toggle" which allows the display to be switched back and forth between seismic and fdom. ',...
        'When both seismic and fdom are shown, there are two clipping controls, the left one being for the ',...
        'seismic and the right one being for the fdom. Feel free to adjust these. Smaller clip ',...
        'numbers mean greater clipping. The word "none" means there is no clipping.  The parameters of the',...
        'computation are shown at right. Hover the pointer over each one for instructions. After changing ',...
        'the parameters push the apply button for a new calculation. There are three different dominant ',...
        'frequency sections, one called Frequency where the amplitudes are the value in Hz of the dominant ',...
        'frequency, another called Bandwidth where the amplitudes are the estimated spectral width in Hz ',...
        'centered on the dominant frequency, and the third called Amplitude where the amplitudes are the ',...
        'spectral amplitude at the dominant frequency. The radio buttons labelled "Display choice" ',...
        'allow you to choose between these. Of these three sections, the Amplitude section is most similar ',...
        'to normal seismic and is the same as evaluating a spectral decomp volume at the dominant frequency ',...
        '(as opposed to displaying a particular frequency in Hz).  However, unlike normal seismic, ',...
        'Amplitude refers to spectral amplitude and hence is never negative. Note that the colorbar ',...
        'always refers to the fdom axes and not the seismic.'];
    hinfo=msgbox(msg,'Instructions for Dominant Frequency tool');
    udat=get(hthisfig,'userdata');
    if(iscell(udat))
        if(isgraphics(udat{1}))
            delete(udat{1});
        end
        udat{1}=hinfo;
    else
        if(isgraphics(udat))
            delete(udat);
        end
        udat=hinfo;
    end
    set(hthisfig,'userdata',udat);
elseif(strcmp(action,'colormap'))
    hcmap=gcbo;
    cmaps=get(hcmap,'string');
    icmap=get(hcmap,'value');
    hcmp=findobj(gcf,'tag','cmapgt');
    switch hcmp.SelectedObject.Tag
        case 'left'
            hax1=findobj(gcf,'tag','seis');
            colormap(hax1,cmaps{icmap});
        case 'right'
            hax2=findobj(gcf,'tag','seisfd');
            colormap(hax2,cmaps{icmap});
        case 'both'
            hax1=findobj(gcf,'tag','seis');
            hax2=findobj(gcf,'tag','seisfd');
            colormap(hax1,cmaps{icmap});
            colormap(hax2,cmaps{icmap});
    end
elseif(strcmp(action,'clip1'))
    hclip=findobj(gcf,'tag','clip1');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    am=udat{2};
    amax=udat{4};
    amin=udat{5};
    sigma=udat{3};
    hax=udat{6}(1);
    if(iclip==1)
        clim=[amin amax];
    else
        clip=clips(iclip);
        clim=[am-clip*sigma,am+clip*sigma];
    end
    set(hax,'clim',clim);
elseif(strcmp(action,'clip2')||strcmp(action,'clip2fromchoice'))
    fromchoice=false;
    if(strcmp(action,'clip2fromchoice'))
        fromchoice=true;
    end
    hresults=findobj(gcf,'tag','results');
    results=get(hresults,'userdata');
    iresult=get(hresults','value');
    choice=nowshowing;
    switch choice
        case 'freq'
            clipdatr=results.clipfd{iresult};
        case 'bw'
            clipdatr=results.clipbwfd{iresult};
        case 'amp'
            clipdatr=results.clipafd{iresult};
    end
    hclip=findobj(gcf,'tag','clip2');
    udat=get(hclip,'userdata');
    if(fromchoice)
        iclip=clipdatr(5);
        set(hclip,'value',iclip');
    else
        iclip=get(hclip,'value');
        clipdatr(5)=iclip;
    end
    clips=udat{1};
    clipdat=udat{2};
    am=clipdat(1);
    amax=clipdat(2);
    sigma=clipdat(4);
    hax=findobj(gcf,'tag','seisfd');
    if(iclip==1)
        %clim=[amin amax];
        clip=max(clips);
    else
        clip=clips(iclip-1);
    end
    %clim=[max([0 am-clip*sigma]),am+clip*sigma];
    clim=[am-clip*sigma,am+clip*sigma];
    set(hax,'clim',clim);
    switch choice
        case 'amp'
            results.clipafd{iresult}=clipdatr;
        case 'freq'
            results.clipfd{iresult}=clipdatr;
        case 'bw'
            results.clipbwfd{iresult}=clipdatr;
    end
    cb=findobj(gcf,'type','colorbar');
    poscb=get(cb,'position');
    delete(cb);
    cb=colorbar;
    cblim=get(cb,'limits');
    if(cblim(1)<0)
        cblim(1)=0;
    end
    if(cblim(2)>amax)
        cblim(2)=amax;
    end
    set(cb,'limits',cblim,'position',poscb)
    set(hresults,'userdata',results);
elseif(strcmp(action,'choice'))
    hresults=findobj(gcf,'tag','results');
    choice=nowshowing;
    iresult=get(hresults,'value');
    results=get(hresults,'userdata');
    hseis2=findobj(gcf,'tag','seisfd');
    hi=findobj(hseis2,'type','image');
    hclip2=findobj(gcf,'tag','clip2');
    switch choice
        case 'amp'
            set(hi,'cdata',results.afd{iresult},'ydata',results.tfd{iresult});
            ud=get(hclip2,'userdata');
            ud{2}=results.clipafd{iresult};
            set(hclip2,'userdata',ud,'value',ud{2}(5));
        case 'freq'

            set(hi,'cdata',results.fd{iresult},'ydata',results.tfd{iresult});
            ud=get(hclip2,'userdata');
            ud{2}=results.clipfd{iresult};
            set(hclip2,'userdata',ud,'value',ud{2}(5));
        case 'bw'
            set(hi,'cdata',results.bwfd{iresult},'ydata',results.tfd{iresult});
            ud=get(hclip2,'userdata');
            ud{2}=results.clipbwfd{iresult};
            set(hclip2,'userdata',ud,'value',ud{2}(5));
    end
    set(hresults,'userdata',results);
    seisplotfdom('clip2fromchoice');  
elseif(strcmp(action,'brighten'))
    hbut=gcbo;
    hbright=findobj(gcf,'tag','brighten');
    if(hbut==hbright)
        inc=.1;
    else
        inc=-.1;
    end
    brighten(inc);
    hbrightness=findobj(gcf,'tag','brightness');
    brightlvl=get(hbrightness,'userdata');
    brightlvl=brightlvl+inc;
    if(abs(brightlvl)<.01)
        brightlvl=0;
    end
    set(hbrightness,'string',['lvl ' num2str(brightlvl)],'userdata',brightlvl)
elseif(strcmp(action,'equalzoom'))
    hbut=gcbo;
    hseis=findobj(gcf,'tag','seis');
    hseissd=findobj(gcf,'tag','seisfd');
    tag=get(hbut,'tag');
    switch tag
        case '1like2'
            xl=get(hseissd,'xlim');
            yl=get(hseissd,'ylim');
            set(hseis,'xlim',xl,'ylim',yl);
            
        case '2like1'
            xl=get(hseis,'xlim');
            yl=get(hseis,'ylim');
            set(hseissd,'xlim',xl,'ylim',yl);
    end
    
elseif(strcmp(action,'close'))
    haveamp=findobj(gcf,'tag','aveamp');
    hspec=get(haveamp,'userdata');
    if(isgraphics(hspec))
        delete(hspec);
    end
    delete(gcf);
elseif(strcmp(action,'hideshow'))
    hbut=gcbo;
    option=get(hbut,'userdata');
    hclip1=findobj(gcf,'tag','clip1');
%     udat1=get(hclip1,'userdata');
    hax1=findobj(gcf,'tag','seis');
    hclip2=findobj(gcf,'tag','clip2');
    %udat2=get(hclip2,'userdata');
    hax2=findobj(gcf,'tag','seisfd');
    hi1=findobj(hax1,'type','image');
    hi2=findobj(hax2,'type','image');
    htoggle=findobj(gcf,'tag','toggle');
    
    switch option
        case 'hide'
            pos1=get(hax1,'position');
            pos2=get(hax2,'position');
            x0=pos1(1);
            y0=pos1(2);
            wid=pos2(1)+pos2(3)-pos1(1);
            ht=pos1(4);
            set(hax1,'visible','off','position',[x0,y0,wid,ht]);
            set(hi1,'visible','off');
            set(hclip1,'visible','off');
            set(hax2,'position',[x0,y0,wid,ht]);
            set(htoggle,'userdata',{pos1 pos2})
            set(hbut,'string','Show seismic','userdata','show')
            set(htoggle,'visible','on');
        case 'show'
            udat=get(htoggle,'userdata');
            pos1=udat{1};
            pos2=udat{2};
            set(hax1,'visible','on','position',pos1);
            set([hi1 hclip1],'visible','on');
            set(hax2,'visible','on','position',pos2);
            set(htoggle','visible','off')
            set(hbut,'string','Hide seismic','userdata','hide');
            set([hi2 hclip2],'visible','on');
    end
elseif(strcmp(action,'toggle'))
    hclip1=findobj(gcf,'tag','clip1');
%     udat1=get(hclip1,'userdata');
    hax1=findobj(gcf,'tag','seis');
    hclip2=findobj(gcf,'tag','clip2');
%     udat2=get(hclip2,'userdata');
    hax2=findobj(gcf,'tag','seisfd');
    hi1=findobj(hax1,'type','image');
    hi2=findobj(hax2,'type','image');
    
    option=get(hax1,'visible');
    switch option
        case 'off'
            %ok, turning on seismic
            set([hax1 hclip1 hi1],'visible','on');
            set([hax2 hclip2 hi2],'visible','off');
            
        case 'on'
            %ok, turning off seismic
            set([hax1 hclip1 hi1],'visible','off');
            set([hax2 hclip2 hi2],'visible','on');
    end
elseif(strcmp(action,'apply'))
    hseis=findobj(gcf,'tag','seis');
    hseis2=findobj(gcf,'tag','seisfd');
    hclip2=findobj(gcf,'tag','clip2');
    hi=findobj(hseis,'type','image');
    seis=get(hi,'cdata');
    t=get(hi,'ydata');
    dt=t(2)-t(1);
    fnyq=.5/(t(2)-t(1));
    hobj=findobj(gcf,'tag','twin');
    val=get(hobj,'string');
    twin=str2double(val);
    if(isnan(twin))
        msgbox('Twin is not recognized as a number','Oh oh ...');
        return;
    end
    if(twin<dt || twin>t(end))
        msgbox(['Twin must be greater than dt and less than ' num2str(t(end))],'Oh oh ...');
        return;
    end
    
    hobj=findobj(gcf,'tag','tinc');
    val=get(hobj,'string');
    tinc=str2double(val);
    if(isnan(tinc))
        msgbox('Tinc is not recognized as a number','Oh oh ...');
        return;
    end
    if(tinc<dt || tinc>twin)
        msgbox('Tinc must be greater than dt and less than Twin','Oh oh ...');
        return;
    end
    hobj=findobj(gcf,'tag','fmax');
    val=get(hobj,'string');
    fmax=str2double(val);
    if(isnan(fmax))
        msgbox('Fmax is not recognized as a number','Oh oh ...');
        return;
    end
    if(fmax<0 || fmax>fnyq)
        msgbox('Fmax must lie between 0 and Nyquist','Oh oh ...');
        return;
    end
    hobj=findobj(gcf,'tag','tfmax');
    val=get(hobj,'string');
    tfmax=str2double(val);
    if(isnan(tfmax))
        msgbox('Tfmax is not recognized as a number','Oh oh ...');
        return;
    end
    if(tfmax<0 || tfmax>fnyq)
        msgbox('Tfmax must lie between t(1) and t(end)','Oh oh ...');
        return;
    end
    t1=clock;
    [fd,afd,bwfd,tfd]=tv_afdom(seis,t,twin,tinc,[fmax tfmax],1,2,1);
    t2=clock;
    timepertrace=round(100000*etime(t2,t1)/size(seis,2))/1000;
    hperf=findobj(gcf,'tag','performance');
    set(hperf,'string',[num2str(timepertrace) ' ms/trace'])
    %create clip info
    iclip=get(hclip2,'value');
    clipfd=[mean(fd(:)),max(fd(:)),min(fd(:)),std(fd(:)) iclip];
    clipafd=[mean(afd(:)),max(afd(:)),min(afd(:)),std(afd(:)) iclip];
    clipbwfd=[mean(bwfd(:)),max(bwfd(:)),min(bwfd(:)),std(bwfd(:)),iclip];
    %determine if amp or frequency or bw
    choice=nowshowing;
    hi=findobj(hseis2,'type','image');
    if(isempty(hi))
        %then we are doing it the first time
        hseis1=findobj(gcf,'tag','seis');
        hi1=findobj(hseis1,'type','image');
        x=get(hi1,'xdata');
        xname=get(get(hseis1,'xlabel'),'string');
        yname=get(get(hseis1,'ylabel'),'string');
        axes(hseis2);
        
        switch choice
            case 'amp'
                imagesc(x,tfd,afd);
                posax=get(hseis2,'position');
                hc=colorbar;
                posc=get(hc,'position');
                set(hseis2,'position',posax);
                set(hc,'position',[posax(1)+posax(3) posc(2:4)]);
                xlabel(xname);ylabel(yname);
                ud=get(hclip2,'userdata');
                ud{2}=clipafd;
                set(hclip2,'userdata',ud);
                set(hseis2,'tag','seisfd');
            case 'freq'
                imagesc(x,tfd,fd);
                posax=get(hseis2,'position');
                hc=colorbar;
                posc=get(hc,'position');
                set(hseis2,'position',posax);
                set(hc,'position',[posax(1)+posax(3)+.1*posc(3) posc(2) .5*posc(3) posc(4)]);
                xlabel(xname);ylabel(yname);
                ud=get(hclip2,'userdata');
                ud{2}=clipfd;
                set(hclip2,'userdata',ud);
                set(hseis2,'tag','seisfd');
            case 'bw'
                imagesc(x,tfd,bwfd);
                posax=get(hseis2,'position');
                hc=colorbar;
                posc=get(hc,'position');
                set(hseis2,'position',posax);
                set(hc,'position',[posax(1)+posax(3)+.1*posc(3) posc(2) .5*posc(3) posc(4)]);
                xlabel(xname);ylabel(yname);
                ud=get(hclip2,'userdata');
                ud{2}=clipbwfd;
                set(hclip2,'userdata',ud);
                set(hseis2,'tag','seisfd');
        end
    else
        switch choice
            case 'amp'
                set(hi,'cdata',afd);
                ud=get(hclip2,'userdata');
                ud{2}=clipafd;
                set(hclip2,'userdata',ud);
            case 'freq'
                set(hi,'cdata',fd);
                ud=get(hclip2,'userdata');
                ud{2}=clipfd;
                set(hclip2,'userdata',ud);
            case 'bw'
                set(hi,'cdata',bwfd);
                ud=get(hclip2,'userdata');
                ud{2}=clipbwfd;
                set(hclip2,'userdata',ud);
        end
    end
    axes(hseis2)
    name=['Fdom for Twin=' num2str(twin) ', Tinc=' num2str(tinc) ', fmax=' num2str(fmax) ', tfmax=' num2str(tfmax)];
    
    %save the results and update hresults
    hresults=findobj(gcf,'tag','results');
    results=get(hresults,'userdata');
    if(isempty(results))
        nresults=1;
        results.names={name};
        results.fd={fd};
        results.afd={afd};
        results.bwfd={bwfd};
        results.tfd={tfd};
        results.twins={twin};
        results.tincs={tinc};
        results.fmaxs={fmax};
        results.tfmaxs={tfmax};
        results.clipfd={clipfd};
        results.clipafd={clipafd};
        results.clipbwfd={clipbwfd};
    else
        nresults=length(results.names)+1;
        results.names{nresults}=name;
        results.fd{nresults}=fd;
        results.afd{nresults}=afd;
        results.bwfd{nresults}=bwfd;
        results.tfd{nresults}=tfd;
        results.twins{nresults}=twin;
        results.tincs{nresults}=tinc;
        results.fmaxs{nresults}=fmax;
        results.tfmaxs{nresults}=tfmax;
        results.clipfd{nresults}=clipfd;
        results.clipafd{nresults}=clipafd;
        results.clipbwfd{nresults}=clipbwfd;
    end
    set(hresults,'string',results.names,'value',nresults,'userdata',results)
    seisplotfdom('clip2');
elseif(strcmp(action,'select'))
    hfig=gcf;
    hresults=findobj(hfig,'tag','results');
    results=get(hresults,'userdata');
    iresult=get(hresults,'value');
    hseis2=findobj(hfig,'tag','seisfd');
    hi=findobj(hseis2,'type','image');
    hamp=findobj(hfig,'tag','amp');
    iamp=get(hamp,'value');
    hclip2=findobj(hfig,'tag','clip2');
    choice=nowshowing;
    switch choice
        case 'amp'
            set(hi,'cdata',results.afd{iresult},'ydata',results.tfd{iresult});
            ud=get(hclip2,'userdata');
            ud{2}=results.clipafd{iresult};
            set(hclip2,'userdata',ud,'value',ud{2}(5));
        case 'freq'
            set(hi,'cdata',results.fd{iresult},'ydata',results.tfd{iresult});
            ud=get(hclip2,'userdata');
            ud{2}=results.clipfd{iresult};
            set(hclip2,'userdata',ud,'value',ud{2}(5));    
        case 'bw'
            set(hi,'cdata',results.bwfd{iresult},'ydata',results.tfd{iresult});
            ud=get(hclip2,'userdata');
            ud{2}=results.clipbwfd{iresult};
            set(hclip2,'userdata',ud,'value',ud{2}(5)); 
    end
    htwin=findobj(hfig,'tag','twin');
    set(htwin,'string',num2str(results.twins{iresult}));
    htinc=findobj(hfig,'tag','twin');
    set(htinc,'string',num2str(results.tincs{iresult}));
    hfmax=findobj(hfig,'tag','fmax');
    set(hfmax,'string',num2str(results.fmaxs{iresult}));
    htfmax=findobj(hfig,'tag','tfmax');
    set(htfmax,'string',num2str(results.tfmaxs{iresult})); 
    seisplotfdom('clip2');
end
end

function [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(data)
% data ... input data
%
% 
% clips ... determined clip levels
% clipstr ... cell array of strings for each clip level for use in popup menu
% clip ... starting clip level
% iclip ... index into clips where clip is found
% sigma ... standard deviation of data
% am ... mean of data
% amax ... max of data
% amin ... min of data

sigma=std(data(:));
am=mean(data(:));
amin=min(data(:));
amax=max(data(:));
nsigma=ceil((amax-amin)/sigma);%number of sigmas that span the data

%clips=linspace(nsigma,1,nclips)';
clips=[20 15 10 8 6 4 3 2 1 .5 .25 .1 .075 .05 .025 .01 .005 .001 .0001]';
if(nsigma<clips(1))
    ind= clips<nsigma;
    clips=[nsigma;clips(ind)];
else
    n=floor(log10(nsigma/clips(1))/log10(2));
    newclips=zeros(n,1);
    newclips(1)=nsigma;
    for k=n:-1:2
        newclips(k)=2^(n+1-k)*clips(1);
    end
    clips=[newclips;clips];
end

clipstr=cell(size(clips));
nclips=length(clips);
clipstr{1}='none';
for k=2:nclips
    clipstr{k}=['clip= ' num2str(sigfig(clips(k),3))];
end
iclip=near(clips,3);
clip=clips(iclip);

end

function choice=nowshowing
hbg=findobj(gcf,'tag','choices');
choice=hbg.SelectedObject.Tag;
end