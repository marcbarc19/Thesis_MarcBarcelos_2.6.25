function [Responses] = SMS_calcBehavInds(ITI,Asyn,IOI,TestRange)

[n,m] = size(Asyn);

for curAsync = 1:m
    
    Abs = abs(Asyn(TestRange,curAsync));

    MeanAsyn(curAsync) = mean(Asyn(TestRange,curAsync));
    StdAsyn(curAsync) = std(Asyn(TestRange,curAsync));
    MABSASYC(curAsync) = mean(Abs);
    STDABSASYC(curAsync) = std(Abs);
    
end

[c,lags] = xcov(IOI(TestRange),ITI(TestRange),1,'coeff'); % cross-correlation analysis for lags -1 0 +1
[d,lags] = xcov(Asyn(TestRange),1,'coeff');

PTI = c(2,:) - c(1,:);
Lag0CC = c(2,:);
Lag1CC = c(1,:);
PTIratio = c(2,:) ./ c(1,:);
PTIweight = PTI / (c(2,:) + c(1,:));
L1ACF_Async = d';

Responses = table(MeanAsyn,StdAsyn,MABSASYC,PTI,Lag0CC,Lag1CC,PTIratio,PTIweight);