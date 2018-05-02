function filtmask=makegaborfilter(t,twin,tinc,gdb,t1,fmin,fmax,fmaxlim,fphase,max_atten)

[tvs,trow,fcol,normf_tout]=fgabor(rand(size(t)),t,twin,tinc,1,gdb,0); %#ok<ASGLU>

filtmask=zeros(length(trow),length(fcol));
for k=1:length(trow)
    if(t1==-1)
        fmn=fmin;
        fmx=fmax;
    else
        fmn=fmin;
        f2=fmax(1)*t1/trow(k);
        if(f2>fmaxlim(1)); f2=fmaxlim(1); end
        if(f2<fmaxlim(2)); f2=fmaxlim(2); end
        fmx=[f2 fmax(2)];
    end
    dt=.5/fcol(end);
    tmax=.5/fcol(2);
    filtmask(k,:)=(filtspec(dt,tmax,fmn,fmx,fphase,max_atten))';
end