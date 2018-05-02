function viewtextheader(trchdr,dname)
% 
% viewtexthdr(trchdr,dname)
%
%

if(nargin<2)
    dname='';
end

ss=get(0,'screensize');

% nlines=40;
fs=10;
% pixperline=10*fs/8;
% fight=nlines*pixperline;%figure height
figwd=1300;
fight=1000;

hfig=figure;
%pos=hfig.Position;
xnot=200;
ynot=ss(4)-fight-200;
hfig.Position=[xnot ynot figwd fight];
uicontrol(hfig,'style','listbox','string',trchdr,'units','normalized','position',[.1 .05 .8 .8],...
    'horizontalalignment','left','fontsize',fs);

if(~isempty(dname))
    set(gcf,'name',['SEGY Text header for ',dname]);
else
    set(gcf,'name','SEGY Text header');
end

