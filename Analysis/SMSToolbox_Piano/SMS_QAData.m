function [TapTime_Int,ITI_Int,asyn_Int,invalid,clean,Interpolated,MissedTap,LargeAsynchrony,NumInterp] = SMS_QAData(NTap,asyn,ITI,IOI,TapTime,BeepTimeRaw,MissedCrit,TestRange,lag)

% Performs quality assurance on SMS data. Checks for invalid trials (missed taps, late (or
% early) taps). Determines whether the trial should be rejected by checking 
% for the presence of consequtive missing taps, or if the number of missed 
% taps is equal to or greater than MissedCrit. Interpolates invalid taps if
% the trial has not been rejected.

% -- INPUT ARGUMENTS --
% NTap - Binary vector of missed taps (1=missed tap)

% asyn - Vector of asynchronies

% ITI - Intertap interval

% IOI - Inter-onset interval

% TapTime - Tap times relative to trial onset

% MissedCrit - Sets number of missed taps to reject trial

% TestTange - indices for the tap/tones of interest



% -- OUTPUT ARGUMENTS --
% TapTime_Int - Interpolated tap times

% ITI_Int - Interpolated intertap intervals

% asyn_Int - Interpolated asynchronies

% invalid - Binary vector of invalid (large asychronies, missed) taps.

% clean - 1 if trial should not be rejected, 0 if trial shold be rejected

% MissedTap - Number of missed taps

% LargeAsynchrony - Number of large Asychnronies

% PropInterp - Proportion of interpolated taps


%% Set output variables
asyn_Int = asyn;
TapTime_Int = TapTime;
ITI_Int = ITI;

%% Determine invalid taps

% Check for large asynchronies
bound = IOI .* .5;
AsyncSc =  (asyn > bound | asyn < -bound);

NTap=NTap==0;
% Check for invalid taps
MissedTap = nnz(NTap(TestRange)); %Number of missed taps within the test range
LargeAsynchrony = nnz(AsyncSc(TestRange)); %Number of large asynchronies within the test range
invalid = NTap + AsyncSc; %Nmuber of invalid trials
NumInterp = sum(invalid(TestRange));

% Report invalid taps
fprintf('-- Checking data quality \n');
fprintf('-- %d Large Asynchronies within Test Range \n',LargeAsynchrony)
fprintf('-- %d Missed Taps within Test Range \n',MissedTap)

%% Provide Diagnostics
% Plot IOI, ITI and Asynchronies.

H1 = figure(1);
%set(H1,'Position',[100,100,800,800])
subplot(2,2,1)
plot(TestRange,ITI(TestRange),TestRange,IOI(TestRange))
title('ITI & IOI')
ylabel('milliseconds')
xlabel('Trial')
legend('ITI','IOI')
subplot(2,2,3)
plot(TestRange,asyn(TestRange),TestRange,bound(TestRange),TestRange,-bound(TestRange))
title('Asynchrony')
ylabel('milliseconds')
xlabel('Trial')

H2 = figure(2);
%set(H2,'Position',[100,1000,1200,200])
scatter(BeepTimeRaw(TestRange),ones(1,length(TestRange)))
hold on
scatter(TapTime(TestRange)-lag,ones(1,length(TestRange)).*.99,[],[1 0 0 ])

axis([0,max([BeepTimeRaw(end),TapTime(end)]),0.9,1.1])

%% Determine whether trial should be rejected
clean = 1; %assume the file is good

%check for too many invalid trials
if sum(invalid) > MissedCrit
    fprintf('-- Warning || Data has %d invalid trials \n',sum(invalid(TestRange)))
    clean = 0; %the file should be rejected
end

%check for consecutive invalid trials
flag = 0;
for curTap = TestRange(1)-1:TestRange(end-1)+1
    if invalid(curTap) + invalid(curTap + 1) == 2
        clean = 0;
        fprintf('-- Warning || Data has consecutive missed taps \n')
        break
        
    end
end

%% Interpolation

Interpolated = 0;
if clean==1 && sum(invalid(TestRange))==0 
fprintf('-- Data is fine \n')
elseif clean==0
fprintf('-- Data REJECTED \n')

elseif clean==1 && sum(invalid(TestRange))>0
    
fprintf('-- Interpolating invalid data using average \n');
Interpolated = 1;
% This function interpolates missing ITI and Asynchronies
% separately

% Find Invalid taps
invalid(1:TestRange(1)-1) = 0;
invalid(TestRange(end)+1:end) = 0;
invInds = find(invalid==1);

% Interpolate missing asynchronies by averageing adjacent asynchronies
theInterps = (asyn(invInds-1) + asyn(invInds+1))./2;
asyn_Int(invInds) = theInterps;

% Interpolate missing TapTimes by averaging adjacent TapTimes
theInterps = (TapTime(invInds+1) - TapTime(invInds-1))./2 + TapTime(invInds-1);
TapTime_Int(invInds) = theInterps;

% Calculate missing ITIs from Interpolated TapTimes
ITI_Int = (TapTime_Int(2:end) - TapTime_Int(1:end-1)).* 1000;
ITI_Int = [0;ITI_Int];

%Report Interpolated Taps
fprintf('-- Check Trials %d \n', invInds);
      

H1 = figure(1);
subplot(2,2,2)
plot(TestRange,ITI_Int(TestRange),TestRange,IOI(TestRange))
title('ITI & IOI')
ylabel('milliseconds')
xlabel('Trial')
legend('ITI','IOI')
subplot(2,2,4)
plot(TestRange,asyn_Int(TestRange),TestRange,bound(TestRange),TestRange,-bound(TestRange))
title('Asynchrony')
ylabel('milliseconds')
xlabel('Trial')

H2 = figure(2);
scatter(TapTime_Int(TestRange)-lag,ones(1,length(TestRange)).*.98,[],[1 0 0 ])

end

disp('Press SpaceBar To Resume Analysis')
pause;
close(H1);
close(H2);
