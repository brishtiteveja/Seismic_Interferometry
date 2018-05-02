function stackg=gabordecon_stackPar(stack,t,twin,tinc,tsmo,fsmo,ihyp,stab,phase,t1,fmin,fmax,fmaxlim,fphase,max_atten)
% GABORDECON_STACKPAR: applies gabor decon to a stacked section, enabled for parallel
%
% stackg=gabordecon_stackPar(stack,t,twin,tinc,tsmo,fsmo,ihyp,stab,phase,t1,fmin,fmax,fmaxlim,fphase,max_atten)
%
% GABORDECON_STACKPAR applies gabordecon to all traces in a stack. The main loop over traces is a
% parfor loop
%
% stack ... stacked section as a matrix of traces. 
% t ... time coordinate for stack
% twin ... half width of gaussian temporal window (sec)
% tinc ... temporal increment between windows (sec)
% tsmo ... size of temporal smoother (sec)
% fsmo ... size of frequency smoother (Hz)
% ihyp ... 1 for hyperbolic smoothing, 0 for ordinary boxcar smoothing
%    Hyperbolic smoothing averages the gabor magnitude spectrum along
%    curves of t*f=constant.
% ************** Default = 1 ***********
% stab ... stability constant
%   ************* Default = 0.000001 **************
% phase ... 0 for zero phase, 1 for minimum phase
%   ************* Default = 1 **************
% 
% stackg ... deconvolved stack
%
% G.F. Margrave, Devon Canada, May 2016
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


if(nargin<9)
    phase=0;
end
if(nargin<8)
    stab=0.000001;
end
if(nargin<7)
    ihyp=1;
end
p=1;
gdb=60;
    
ntr=size(stack,2);

if(length(t)~=size(stack,1))
    error('invalid t coordinate vector')
end

stackg=zeros(size(stack));

small=100*eps;

% t0=clock;
filtmask=makegaborfilter(t,twin,tinc,gdb,t1,fmin,fmax,fmaxlim,fphase,max_atten);
parfor k=1:ntr
    tmp=stack(:,k);
    if(sum(abs(tmp))>small)%avoid deconvolving a zero trace
        stackg(:,k)=gabordeconfilt(stack(:,k),t,twin,tinc,tsmo,fsmo,ihyp,stab,phase,p,gdb,filtmask);
    end
end

% tnow=clock;
% time_used=etime(tnow,t0);
% time_per_trace=time_used/ntr;
% disp(['Finished Gabor, total time= ' num2str(time_used/60) ' min'])
% disp(['time-per-trace= ' int2str(1000*time_per_trace) ' ms'])

amax=max(abs(stackg(:)));
stackg=stackg/amax;
    

