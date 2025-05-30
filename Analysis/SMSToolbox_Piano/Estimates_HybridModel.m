function [m,alpha,sT,sM, LL]=Estimates_HybridModel(s,r,e)
%%%% Function for estimating parameters of the hybrid model

%           s = IOI
%           r = ITI
%           e = asyn

%%% output:
%           alphaE = phase correction parameter
%           mE = prediction/tracking parameter
%           st = timekeeper noise estimate
%           sm = motor noise estimate 
%           LL = log likelihood -- fit of the model to the data


% b-Ax : check simulations & new formulae
b=r(5:end)-s(4:(end-1)); % vector
A=nan(length(b),2); % matrix
A(:,1)=s(4:(end-1))-s(3:(end-2))-e(4:(end-1));
A(:,2)=-e(4:(end-1));

% b=r(4:end)-s(3:(end-1)); % vector
% A=nan(length(b),2); % matrix
% A(:,1)=s(3:(end-1))-s(2:(end-2))-e(3:(end-1));
% A(:,2)=-e(3:(end-1));
%ddd=b-A*[m;alpha]-noise(2:(end-1)); % check simulations & new formulae
%figure(1);clf;plot(ddd);

% correct mean for bGLS
b=b-mean(b);
A(:,1)=A(:,1)-mean(A(:,1));
A(:,2)=A(:,2)-mean(A(:,2));

% do bGLS
[x,sM,sT,LL]=doGLS2bound_2014simplified(A,b);

% calculate estimates parameters
m=x(1);
alpha=x(2);
