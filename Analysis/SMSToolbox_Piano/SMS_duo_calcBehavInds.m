function [Responses] = SMS_calcBehavInds(ITI,Asyn,IOI,TestRange)

[n,m] = size(Asyn);

for curAsync = 1:m
    
    Abs = abs(Asyn(TestRange,curAsync));

    MeanSigAsyn(curAsync) = mean(Asyn(TestRange,curAsync));
    StdSigAsyn(curAsync) = std(Asyn(TestRange,curAsync));
    MeanAbsAsyn(curAsync) = mean(Abs);
    StdAbsAsyn(curAsync) = std(Abs);
    
end

[CC,lags] = xcov(ITI(TestRange,1),ITI(TestRange,2),1,'coeff'); % cross-correlation analysis for lags -1 0 +1
%[d,lags] = xcov(Asyn(TestRange),1,'coeff');

% PTI = c(2,:) - c(1,:);
% Lag0CC = c(2,:);
% Lag1CC = c(1,:);
% PTIratio = c(2,:) ./ c(1,:);
% PTIweight = PTI / (c(2,:) + c(1,:));
% L1ACF_Async = d;
% 

CC=CC';

Responses = table(MeanSigAsyn,StdSigAsyn,MeanAbsAsyn,StdAbsAsyn,CC);