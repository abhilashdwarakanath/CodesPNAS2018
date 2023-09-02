function Pout=jitter2(Pin,isirand,dt);

% Return a 'jittered' version of the process
% Pin.  isirand is a function handle that returns
% random number from the jitter distribution.
% dt is the binsize of the (to determine units).
% Pout=jitter(Pin,jrand,dt)
% To impose a Gaussian (normally distributed) 
% jitter, you might call:
%      pout=jitter(pin,@randn,dt);
% This script is optimised to run with parfor. 
% dt is the amount the data needs to be jittered by in number of samples.
% Abhilash Dwarakanath. MPI biological cybernetics

% Initialize Pout
n=numel(Pin);
Pout=zeros(1,n);

dt=dt/1000;

% Get a representation of Pin 
% in terms of spike times
s=find(Pin);


% Fill Pout making sure not to
% use indices that are out of bounds
for i=1:numel(s)
  
    r = randn;
    %J=round(s(i)+(feval(isirand)+dt));
	J=round(s(i)+(feval(isirand)/dt));
    if and(J>0, J<n)==1
    Pout(J)=Pout(J)+1;
%     else
%         Pout(1) = Pout(1)+1;
%     end
%     if J<n
%         Pout(J)=Pout(J)+1;
%     else
%         Pout(end)=Pout(end)+1;
%     end
    end
   
end

Pout(Pout==2)=1;
    


