function [gammaE,mE,alphaE,stE,smE,LLE]=Estimates_JointModelAlpha(s,r,e,L,H)
%%%% Function for estimating parameters of the Joint Model Alpha.

%%%% Parameterspace for alpha is restricted in order to obtain better estimates. 
%%%% Motor noise is set to zero [using doGLS2bound_2014simplified_joint_nomotor]
%%%% to allow comparison between models using LL. Models therefore need to have the
%%%% same number of parameters.In order to obtain motor noise parameters
%%%% comment line 35 [sM=0; %NOOOTE I KILLED HERE MOTOR VARIANCE] in
%%%% function doGLS2bound_2014simplified_joint_nomotor

%%% input:  
%           s = IOI
%           r = ITI
%           e = asyn
%           L = lower bound of the restriction (we used -0.8)
%           H = upper bound of the restriction (we used -0.1)

%%% output:
%           alphaE = phase correction parameter
%           gammaE = anticipatory error correction parameter
%           mE = prediction/tracking parameter
%           st = timekeeper noise estimate
%           sm = motor noise estimate --> set to zero for comparison using
%           the LL.
%           LL = log likelihood -- fit of the model to the data


% M.C. van der Steen & T.A.B. Dollevoet & N. Jacoby

% b-Ax
b=r(5:end);
A(:,1)=s(4:(end-1))-s(3:(end-2));
A(:,2)=s(4:(end-1))-e(4:(end-1));
A(:,3)=-e(4:(end-1));

% correct for mean bGLS
bm=b-mean(b);
Am=A;
Am(:,1)=A(:,1)-mean(A(:,1));
Am(:,2)=A(:,2)-mean(A(:,2));
Am(:,3)=A(:,3)-mean(A(:,3));

% do bGLS L = lower bound alpha / H = upper bound alpha
[xB,sMB,sTB,LLE]=doGLS2bound_2014simplified_joint_nomotor(Am,bm,L,H);

% calculate parameter estimates
gammaE=1-xB(2);
mE=xB(1)/(1-gammaE);
alphaE=xB(3)/gammaE;
stE=sTB/sqrt(1+2*(gammaE^2)-2*gammaE);
smE=sMB;