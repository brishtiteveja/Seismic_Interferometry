function viewtraceheaders(tracehdrs,msg)
% VIEWTRACEHEADERS - Put up a GUI to browse SEGY trace headers
%
% viewtraceheaders(tracehdrs,msg)
%
% tracehdrs ... trace header structure as returned from SegyTrace or SegyFile or readsegy.
% msg ... text message to display at the top of the window
% ************ default '' *************
%
% G.F. Margrave, Devon, 2017
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
% 
if(~ischar(tracehdrs))
    action='init';
else
    action=tracehdrs;
end

if(strcmp(action,'init'))
    
    if(nargin<2)
        msg='';
    end
    
    fn=fieldnames(tracehdrs);
    v1=tracehdrs.(fn{1});
    
    ntraces=length(v1);
    nfields=length(fn);
    
    %determine segyrev
    segyrev=0;
    if(isfield(tracehdrs,'CdpX'))
        segyrev=1;
    end
    
    data=cell(nfields,4);
    ktrace=1;
    for k=1:nfields
        byte=traceword2byte(fn{k},segyrev);
        data{k,4}=byte;
        data{k,2}=class(tracehdrs.(fn{k})(ktrace));
        data{k,3}=1;%just a placeholder
        data{k,1}=fn{k};
    end
    
    
    
    hf=figure;
    pos=get(hf,'position');
    figwid=600;
    fight=800;
    
    set(hf,'name','Trace Header Viewer','numbertitle','off','menubar','none','toolbar','none',...
        'position',[pos(1), 200, figwid fight]);
    %'position',[pos(1)+.5*(pos(3)-figwid), pos(2)+.5*(pos(4)-fight), figwid fight]);
    
    xnot=.1;
    ynot=.9;
   
    wid=.1;ht=.02;
    xsep=.05;
    ysep=.02;
    xnow=xnot;
     ynow=ynot+2*ht;
    uicontrol(hf,'style','text','string',msg,'tag','msg','units','normalized',...
        'position',[xnow,ynow,.8,ht],'fontsize',10,'fontweight','bold');
    ynow=ynot;
    
    uicontrol(hf,'style','text','string','Trace number','units','normalized','position',...
        [xnow ynow wid ht]);
    uicontrol(hf,'style','edit','string',int2str(ktrace),'tag','tracenum','units','normalized',...
        'position',[xnow+wid, ynow, wid ht],'userdata',{tracehdrs, segyrev, fn, ntraces, nfields, data},...
        'callback','viewtraceheaders(''newtrace'');');
    xnow=xnow+2*wid+xsep;
    uicontrol(hf,'style','pushbutton','string','Next','tag','next','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','viewtraceheaders(''newtrace'');');
    xnow=xnow+wid;
    uicontrol(hf,'style','pushbutton','string','Prev','tag','prev','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','viewtraceheaders(''newtrace'');');
    xnow=xnow+wid+xsep;
    uicontrol(hf,'style','text','string','Incr:','units','normalized',...
        'position',[xnow,ynow,wid,ht]);
    xnow=xnow+wid;
    uicontrol(hf,'style','edit','string','1','tag','inc','units','normalized',...
        'position',[xnow,ynow,wid,ht]);
    xnow=xnot;
    ynow=ynow-ht-ysep;
    uicontrol(hf,'style','text','string',['Total traces: ' int2str(ntraces)],'units','normalized','position',...
        [xnow ynow 3*wid ht]);
    xnow=xnot+3*wid;
    uicontrol(hf,'style','text','string',['SEGY revision #' int2str(segyrev)],'units','normalized','position',...
        [xnow ynow 3*wid ht]);
    
    panht=.8;panwid=.8;
    xnow=xnot;ynow=ynow-panht-ysep;
    hdp=uipanel(hf,'tag','data_panel','units','normalized','position',[xnow,ynow,panwid,panht]);
    
    tabht=3;
    tabwid=1;
    yn=-2;
    xn=0;
    fudge=.96;
    colwids={.3*figwid*panwid*fudge, .3*figwid*panwid*fudge, .2*figwid*panwid*fudge, .2*figwid*panwid*fudge};
    
    data=filldata(ktrace);
    
    ht=uitable(hdp,'data',data,'ColumnName',{'Name','Format','Value','Byte'},'units','normalized',...
        'position',[xn yn tabwid tabht],'columnwidth',colwids,'rowname',[]);
    
    uicontrol(hf,'style','slider','units','normalized','tag','slider',...
        'position',[xnow+panwid,ynow,.25*wid,panht],'value',1,'callback',{@sliderstep,ht});
    
    return;
elseif(strcmp(action,'newtrace'))
    hbut=gcbo;
    hnext=findobj(gcf,'tag','next');
    hprev=findobj(gcf,'tag','prev');
    htrace=findobj(gcf,'tag','tracenum');
    udat=get(htrace,'userdata');
    ntraces=udat{4};
    hinc=findobj(gcf,'tag','inc');
    ht=findobj(gcf,'type','uitable');
    if(hbut==hnext)
        sgn=1;
    elseif(hbut==hprev)
        sgn=-1;
    else
        sgn=0;
    end
    inc=str2double(get(hinc,'string'));
    if(isnan(inc))
        inc=1;
        set(hinc,'string','1');
    end
    
    ktrace=str2double(get(htrace,'string'));
    if(isnan(ktrace))
        ktrace=1;
        set(htrace,'string','1');
    end
    
    ktrace=ktrace+sgn*inc;
    if(ktrace>ntraces)
        ktrace=ntraces;
    end
    if(ktrace<1)
        ktrace=1;
    end
    
    set(htrace,'string',int2str(ktrace));
    
    ht.Data=filldata(ktrace); 
    
end



function data=filldata(ktrace)

htrace=findobj(gcf,'tag','tracenum');
udat=get(htrace,'userdata');
tracehdrs=udat{1};
fn=udat{3};
nfields=udat{5};
data=udat{6};


for k=1:nfields
    data{k,3}=tracehdrs.(fn{k})(ktrace);
end

function sliderstep(src,eventdata,arg1) %#ok<INUSL>
val=get(src,'value');
set(arg1,'position',[0 -2*val 1 3]);