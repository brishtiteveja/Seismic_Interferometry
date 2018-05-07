function vec=strcell2num(scell)
% STRCELL2NUM: convert a cell array of strings to a vector of numbers
%
% scell = num2strcell(vec,sigfig)
%
% vec ... vector of numbers
% sigfig ... number of significant figures desired. Is -1 then the vec is considered to be times
% which must be displayed to the nearest millisecond.
% ************** default =4 **************
%
% Converts a vector of numbers to a cell vector of strings. Thus vec(k)=str2double(scell{k}).
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

if(nargin<2)
    sigfig=4;
end

scell=cell(size(vec));
for k=1:length(vec)
    if(sigfig==-1)
        scell{k}=time2str(vec(k));
    else
        scell{k}=num2str(vec(k),sigfig);
    end
end