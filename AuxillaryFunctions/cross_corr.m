function [corrfunc,lags]=cross_corr(p1,p2,winrad,dt)
% Compute the cross correlation function of processes p1 and p2
% dt is the binsize used in the representation of p1 and p2
% winsize is the window radius over which to compute the CCG
% H is the cross correlation function and X is the time domain
% of the function (X=-winrad:dt:winrad).
% [corrfunc,lags]=cross_corr(p1,p2,winrad,dt)

n1=numel(p1);
n2=numel(p2);
%if n1~=n2
%error('Inputs have different sizes.');%
%end

%n=n1;    % n is the number of time units in each process (the size of the vectors)
%T=n*dt;  % T is the length of the spike train in the correct units


% Change the units of winrad to the binsize
winrad=floor(winrad/dt);


% Estimate the rates of each process
%n1=sum(p1)/dt;
%n2=sum(p2)/dt;


% Initialize the vector that will hold the correlation function
corrfunc=zeros(1,2*winrad+1);

% Get the spike times of p2
if length(p1) > length(p2)
    s=find(p2);
    % Fill the vector un-normalized correlation function vector
    p1_temp=[zeros(1,winrad) p1 zeros(1,winrad)];      % stick some zeros on either end of p1 for easier indexing
    for i=1:numel(s)
        corrfunc=corrfunc+p1_temp(s(i):s(i)+2*winrad);   %p1_temp(ti:ti+2*winrad) is really just p1(ti-winrad:ti+winrad) with zeros for the non-existent entries
    end
    
else
    s=find(p1);
    % Fill the vector un-normalized correlation function vector
    p2_temp=[zeros(1,winrad) p2 zeros(1,winrad)];      % stick some zeros on either end of p1 for easier indexing
    for i=1:numel(s)
        corrfunc=corrfunc+p2_temp(s(i):s(i)+2*winrad);   %p1_temp(ti:ti+2*winrad) is really just p1(ti-winrad:ti+winrad) with zeros for the non-existent entries
    end
end



% Normalise to get the correlation function. This is the standard normalisation we use. Depending on your needs, change the value of normalisation.
%corrfunc=(corrfunc-n1*n2/n)./(T*dt);


% Return lags
lags=-winrad*dt:dt:winrad*dt;
