function [gammaE,mE,betaE,stE,smE,LLE]=Estimates_JointModelBeta(s,r,e,L,H)
%%%% Function for estimating parameters of the Joint Model Beta.

%%%% Parameterspace for beta is restricted in order to obtain better estimates. 
%%%% Motor noise is set to zero [using doGLS2bound_2014simplified_jointbeta_nomotornoise]
%%%% to allow comparison between models using LL. Models therefore need to have the
%%%% same number of parameters. In order to obtain motor noise parameters
%%%% comment line 34 [sM=0; %NOOOTE I KILLED HERE MOTOR VARIANCE] in
%%%% function doGLS2bound_2014simplified_jointbeta_nomotornoise

%%% input:  
%           s = IOI
%           r = ITI
%           e = asyn
%           L = lower bound of the restriction (we used 0)
%           H = upper bound of the restriction (we used 1.1)

%%% output:
%           betaE = period correction parameter
%           gammaE = anticipatory error correction parameter
%           mE = prediction/tracking parameter
%           st = timekeeper noise estimate
%           sm = motor noise estimate --> set to zero for comparison using
%           the LL.
%           LL = log likelihood -- fit of the model to the data


% M.C. van der Steen & T.A.B. Dollevoet & N. Jacoby

e=e-mean(e);
% b-Ax
b=r(5:end);
A(:,1)=s(4:(end-1))-s(3:(end-2));
A(:,2)=s(4:(end-1))-e(4:(end-1));
A(:,3)=-cumsum(e(4:(end-1)));

% correct for mean bGLS
bm=b-mean(b);
Am=A;
Am(:,1)=A(:,1)-mean(A(:,1));
Am(:,2)=A(:,2)-mean(A(:,2));
% Am(:,3)=A(:,3)-mean(A(:,3));

% do bGLS L = lower bound alpha / H = upper bound alpha
[xB,sMB,sTB,LLE]=doGLS2bound_2014simplified_jointbeta_nomotornoise(Am,bm,L,H);

% calculate parameter estimates
gammaE=1-xB(2);
mE=xB(1)/(1-gammaE);
betaE=xB(3)/gammaE;
stE=sTB/sqrt(1+2*(gammaE^2)-2*gammaE);
smE=sMB;
    