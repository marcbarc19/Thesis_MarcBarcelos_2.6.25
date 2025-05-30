%% Thesis ADAM Analysis
% February 2025

% This script analyzes jazz trio performance data provided by Cheston et al., 2024.
% The script organises and stores data and the parameters used to run the analysis. 
% It first opens the target file, performs quality assurance, interpolates missing data,
% calculates behavioural indicies, and then estimates parameters from
% the Adaptation and Anticipation Model (ADAM) from Keller & Appel (2010) for each file.
%
% ------
% NB: This script calls functions and scripts from 'SMSToolbox_Piano', which should be findable in
% Matlab path. The script was created and tested in MATLAB R2022a.
% ------

% ------
% INPUT
% 1. Input data should be in .csv files with 5 columns (ITI, TapTime, ASYN, IOI, ToneTime):
%
% ... ITI = inter-tone intervals for TARGET musician (the one for which
% parameters are being estimated)
%
% ... ToneTime = Keystroke time for TARGET pianist (the one for which
% parameters are being estimated)
%
% ... ASYN = time difference between target pianist keystroke and reference
% pianist keystroke (target - reference; negative asynchronies indicate
% target keystrokes earlier than reference)
%
% ... IRI = inter-response intervals for REFERENCE musician (the one for which
% parameters are NOT being estimated)
%
% ... ResponseTime = Keystroke time for REFERENCE pianist (the one for which
% parameters are NOT being estimated)

% ------
% HOW TO USE 
% 1. Set path to folder where all data files are located at 'cd' (if
% changing current directory not wanted, use 'datadir' option)
%
% 2. Set the row number of the data file where analyses should start
% ('startRow') and end ('endRow')
%
% 3. Specify trial rejection criteria: 
% - Total number of missing keystrokes allowed ('MissedCrit'); 
% - Total number of consective missing keystrokes allowed; 
% - Criterion for deleting outlier asynchronies
%
% 4. Set ADAM model parameters

% ------
% OUTPUT
% Analysis results are output in .csv file in same folder as data files
% Contents of the output file are described at the end of the script.


%% INPUT ANALYSIS PARAMETERS
clc

% FILE INFORMATION
cd 'C:\Users\marc1\OneDrive\Desktop\School\Aarhus University\Current Classes\Thesis\Data\processed_data\concat_data' % Path to directy with all data files
datadir = cd;

startRow = 1; % How many tones to omit at start (NB: This will be the first tone to compute ITI, so final data series will start at startRow+1)
endRow = 0; % How many tones to omit at end

% REJECTION CRITERIA
MissedCrit = 250; % Permitted number of total missing events (taps & outlier asynchronies)
ConsecMiss = 1; % Permitted number of consecutive missing events (taps & outlier asynchronies)
%TestRange = 6:endRow-1; % include only taps within this range for quality assurance and analysis
CritAsync = 3; % asynchronies greater than this will be deemed invalid (in beats)

% MODEL PARAMETERS
% Mod_method = {'adapt';'jointBeta'}; % Versions of ADAM model to run
ITER = 20;
Kratio = 1;
LowBound = -1.1;   % Lower Bound, for some parameters 
UpBound = 2.1;   % Upper bound, for some parameters
% LowBound = -1;   % Lower Bound, for some parameters (PK changed for test purposes)
% UpBound = 2;   % Upper bound, for some parameters (PK changed for test purposes)

% STORE PARAMETERS

Parameters.startRow = startRow;
Parameters.endRow = endRow;
Parameters.MissedCrit = MissedCrit;
% Parameters.TestRange = TestRange;
Parameters.CritFactor = CritAsync;
% Parameters.Mod_method = Mod_method;
Parameters.ITER = ITER;
Parameters.Kratio = Kratio;
Parameters.LowBound = LowBound; 
Parameters.UpBound = UpBound;   

SMS_params.Parameters = Parameters;

%% PREPARE DATA
% Get files
FileList = dir(fullfile(datadir, '*.csv')); % Make list of files in data directory

% Split file names so that participant number, trial number, and condition
% can be extracted
metadataDir = 'C:\Users\marc1\OneDrive\Desktop\School\Aarhus University\Current Classes\Thesis\Data\processed_data\metadata'; % Define metadata directory

for FileN = 1:length(FileList)
% Extract filename and split components
SplitList(FileN,:) = strsplit(FileList(FileN,1).name, {'_','.'}, 'CollapseDelimiters', true);

% Extracting relevant parts for matching
curID(FileN,:) = SplitList(FileN,1);
curSong(FileN,:) = SplitList(FileN,2);
curArtist(FileN,:) = SplitList(FileN,3);

% Create base filename (including only song and artist)
baseFilename = strjoin({curSong{FileN}, curArtist{FileN}}, '_'); % Join curSong and curArtist

% Find the metadata file (ignore numeric prefixes)
metadataFiles = dir(fullfile(metadataDir, [baseFilename, '_metadata.json']));

if ~isempty(metadataFiles) % If at least one matching metadata file exists
    metadataFile = fullfile(metadataDir, metadataFiles(1).name); % Take the first match
    
    % Read JSON file and extract tempo
    jsonData = fileread(metadataFile);
    jsonStruct = jsondecode(jsonData);
    if isfield(jsonStruct, 'tempo')
        tempo(FileN,1) = jsonStruct.tempo; % Store tempo value
    else
        tempo(FileN,1) = NaN; % Store NaN if "tempo" field is missing
    end
else
    tempo(FileN,1) = NaN; % Store NaN if metadata file is missing
end


FileName = fullfile(datadir,FileList(FileN).name);
rawData = readtable(FileName);
testData = rawData(startRow:height(rawData)-endRow,:); 

ITIOrig = table2array(testData(:,2)); % Inter-keystroke interval for target pianist (being estimated); may be interpolated
TapOrig = table2array(testData(:,1)); % Keystroke data for target pianist (being estimated); may be interpolated
AsyncOrig = table2array(testData(:,3)); % Asynchrony data from Max; will be checked for outliers and may be interpolated
IOIOrig = table2array(testData(:,5)); % IOI data for reference pianist (not being estimated)
ToneOrig = table2array(testData(:,4)); % Tone onset data for reference pianist (not being estimated)

%% DATA CLEANING

% Check for missing taps (keystroke data for target pianist; the one being estimated)
TapClean = TapOrig;
%TapClean(TapClean==0) = NaN; % Codes missing taps as NaN
nanTap = isnan(TapClean); % NaN coded as 1
nNaNTap(FileN,:) = sum(nanTap(2:end,:)); %!!! Counts total number of missing taps before interpolation (just for portion of trial that will be analyzed)

% Check for missing tones (keystroke data for reference pianist; the one NOT being estimated)
ToneClean = ToneOrig;
nanTone = isnan(ToneClean); % NaN coded as 1
nNaNTone(FileN,:) = sum(nanTone(2:end,:)); %!!! Counts total number of missing tones before interpolation (just for portion of trial that will be analyzed)

% Check for asynchrony outliers
AsyncOut = AsyncOrig;
AsyncOut(abs(AsyncOut)>CritAsync) = NaN; % Codes taps where asynchrony > citical value as NaN
nanAsync = isnan(AsyncOut); % NaN coded as 1
nNaNAsync(FileN,:) = sum(nanAsync(2:end,:)); %!!! Counts total number of missing & outlier asynchronies
nAsyncOutlier(FileN,:) = nNaNAsync(FileN,:) - nNaNTap(FileN,:) - nNaNTone(FileN,:); %!!! Computes total number of missing & outlier asynchronies (by subtracting N missing Taps & Tones)

% Counts total number of missing events based on asynchronies
nanAllCombined = nanAsync;
nanCombined = nanAllCombined > 0; % Does logical test needed for interpolation
nNaNCombined(FileN,:) = sum(nanCombined); %!!! Counts total number of missing events before interpolation (just for portion of trial that will be analyzed)

consecNans(FileN,:) = max(nanCombined(1:end-1)+nanCombined(2:end)) > ConsecMiss+1; %!!! Tests whether there are any consecutive missing taps 
TooManyMissing(FileN,:) = nNaNCombined(FileN,:) > MissedCrit; %!!! Tests whether number of missing taps > critical value setting

% Interpolate if not too many missing
if TooManyMissing(FileN,:) == 0 && consecNans(FileN,:) == 0
    % ------Tap Interpolation------
    Tap = TapClean;  % Interpolates if necessary
    tTap = 1:numel(Tap);
    tTap = tTap';

    nanTap = isnan(TapClean);
    Tap(nanTap) = interp1(tTap(~nanTap), TapClean(~nanTap), tTap(nanTap), 'linear');
    
    % If the first value is NaN
    if isnan(Tap(1))  % If the first value is NaN
        firstValidIndexTap = find(~nanTap, 1, 'first');  % Find the first valid value
        if ~isempty(firstValidIndexTap)
            Tap(1:firstValidIndexTap-1) = Tap(firstValidIndexTap);  % Fill NaNs before the first valid value
        end
    end

    % If the last value is NaN
    if isnan(Tap(end))  % If the last value is NaN
        lastValidIndexTap = find(~nanTap, 1, 'last');  % Find the last valid value in Tap
        if ~isempty(lastValidIndexTap)
            endIndexTap = min(lastValidIndexTap + ConsecMiss, numel(Tap));  % Determine the range of NaNs to be filled with previous value, limited by ConsecMiss
            Tap(lastValidIndexTap+1:endIndexTap) = Tap(lastValidIndexTap);  % Fill the permitted number of NaNs
        end
    end

    % ------Tone Interpolation------
    Tone = ToneClean;  % Interpolates if necessary
    tTone = 1:numel(Tone);
    tTone = tTone';

    nanTone = isnan(ToneClean);
    Tone(nanTone) = interp1(tTone(~nanTone), ToneClean(~nanTone), tTone(nanTone), 'linear');
    
    % If the first value is NaN
    if isnan(Tone(1))  % If the first value is NaN
        firstValidIndexTone = find(~nanTone, 1, 'first');  % Find the first valid value
        if ~isempty(firstValidIndexTone)
            Tone(1:firstValidIndexTone-1) = Tone(firstValidIndexTone);  % Fill NaNs before the first valid value
        end
    end

    % If the last value is NaN
    if isnan(Tone(end))  % If the last value is NaN
        lastValidIndexTone = find(~nanTone, 1, 'last');  % Find the last valid value in Tone
        if ~isempty(lastValidIndexTone)
            endIndexTone = min(lastValidIndexTone + ConsecMiss, numel(Tone));  % Determine the range of NaNs to be filled with previous value, limited by ConsecMiss
            Tone(lastValidIndexTone+1:endIndexTone) = Tone(lastValidIndexTone);  % Fill the permitted number of NaNs
        end
    end



% ------Compute new asynchronies & ITIs (and trim IOI to match size)------
Async = Tap(2:end-1) - Tone(2:end-1); % Don't include first and last event in case interpolation failed for these
ITI = Tap(2:end-1)-Tap(1:end-2);
IOI = Tone(2:end-1)-Tone(1:end-2);

%% DESCRIPTIVE STATS
N_events(FileN,:) = length(Async);

% Compute behavioral performance measures
Mn_Signed_Async(FileN,:) = mean(Async);
Mn_Abs_Async(FileN,:) = mean(abs(Async));
Mn_ITI(FileN,:) = mean(ITI);
SD_Signed_Async(FileN,:) = std(Async);
SD_Abs_Async(FileN,:) = std(abs(Async));
SD_ITI(FileN,:) = std(ITI);
CV_Signed_Async(FileN,:) = std(Async)/mean(ITI);
CV_Abs_Async(FileN,:) = std(abs(Async))/mean(ITI);
CV_ITI(FileN,:) = std(ITI)/mean(ITI);

if ~isempty(Async)
    % Calculate the smallest and largest differences using Tap and Tone
    Min_Tap_Tone = min([Tap(1), Tone(1)]); % Smallest Tap or Tone value at the first event
    Max_Tap_Tone = max([Tap(end), Tone(end)]); % Largest Tap or Tone value at the last event
    
    % Normalize using the difference between max and min, divided by N_events
    Normalized_Value = (Max_Tap_Tone - Min_Tap_Tone) / N_events(FileN,:);
    
    % Compute the new descriptive statistic
    Norm_SD_Signed_Async(FileN,:) = (SD_Signed_Async(FileN,:) / Normalized_Value) * 10;
else
    % Handle cases where Async is empty
    Norm_SD_Signed_Async(FileN,:) = NaN;
end

% Calculate the frequency of each unique IOI
[IOI_values, ~, idx] = unique(IOI);
IOI_counts = histc(idx, 1:numel(IOI_values));

% Normalize to get probabilities
IOI_probs = IOI_counts / sum(IOI_counts);

% Compute Shannon entropy
entropy_interval = -sum(IOI_probs .* log2(IOI_probs));
INT_ent(FileN,:) = entropy_interval;

% Calculate transitions (pairwise IOIs)
transitions = [IOI(1:end-1), IOI(2:end)];

% Find unique pairs and calculate their counts
[transition_values, ~, transition_idx] = unique(transitions, 'rows', 'stable');
transition_counts = histc(transition_idx, 1:size(transition_values, 1));

% Normalize to get probabilities
transition_probs = transition_counts / sum(transition_counts);

% Compute Shannon entropy for transitions
entropy_IOI = -sum(transition_probs .* log2(transition_probs));
IOI_ent(FileN,:) = entropy_IOI;


else
N_events(FileN,:) = length(Async);

Mn_Signed_Async(FileN,:) = NaN;
Mn_Abs_Async(FileN,:) = NaN;
Mn_ITI(FileN,:) = NaN;
SD_Signed_Async(FileN,:) = NaN;
SD_Abs_Async(FileN,:) = NaN;
SD_ITI(FileN,:) = NaN;
CV_Signed_Async(FileN,:) = NaN;
CV_Abs_Async(FileN,:) = NaN;
CV_ITI(FileN,:) = NaN;
IOI_ent(FileN,:) = NaN;
INT_ent(FileN,:) = NaN;

end

%% ADAM
% % Compute ADAM parameter estimates
% % 
% % Note that JointModelBeta assumes '1 - gamma' is anticipatory error
% % correction weight, so higher gamma leads to lower anticipatory error
% % correction. This can be 'corrected' in analyses by first transforming
% % estimates using: 1/gamma
% 
if TooManyMissing(FileN,:) == 0 && consecNans(FileN,:) == 0

s = IOI;
r = ITI;
e = Async;

TestRange = 1:length(e);

% % Adaptation Only version [alphaE,betaE,sTE,sME,LL]
AdaptParams(FileN,:) = SMS_fitModel(SMS_params.Parameters, r, e, s, 'adapt', TestRange);

% % Full Joint version [gammaE,mE,betaE,sTE,sME,LLE]
JointParams(FileN,:) = SMS_fitModel(SMS_params.Parameters, r, e, s, 'jointBeta', TestRange);

else
    AdaptParams(FileN,:) = table(NaN, NaN, NaN, NaN, NaN, NaN, NaN);
    JointParams(FileN,:) = table(NaN, NaN, NaN, NaN, NaN, NaN, NaN);
end

% Save critical values used in pre-processing
MissedCriterion(FileN,:) = MissedCrit;
ConsecutiveMiss(FileN,:) = ConsecMiss;
CriterionAsync(FileN,:) = CritAsync;

end

%% COMPILE & SAVE OUTPUT
    % Estimation results
    Results.File = table(curID, curSong, curArtist);
    Results.File.Properties.VariableNames = ["ID","Song","Artist"];
    Results.BehavMeasures = table(Mn_Signed_Async, Mn_Abs_Async, Mn_ITI, Norm_SD_Signed_Async, SD_Signed_Async, SD_Abs_Async, SD_ITI, CV_Signed_Async, CV_Abs_Async, CV_ITI);
    Results.AdaptParams = AdaptParams;
    Results.AdaptParams.Properties.VariableNames = ["Adapt_Alpha","Adapt_Beta","Adapt_Delta","Adapt_Gamma","Adapt_Tn","Adapt_Mn","Adapt_LLE"];
    Results.JointParams = JointParams;
    Results.JointParams.Properties.VariableNames = ["Joint_Alpha","Joint_Beta","Joint_Delta","Joint_Gamma","Joint_Tn","Joint_Mn","Joint_LLE"];
    Results.Entropy = table(INT_ent, IOI_ent);
    Results.Entropy.Properties.VariableNames = ["Interval Entropy", "IOI Entropy"];
    Results.Missing = table(tempo, TooManyMissing, consecNans, nAsyncOutlier, nNaNTap, nNaNTone, nNaNCombined, N_events, MissedCriterion, CriterionAsync, ConsecutiveMiss);
    Results.Missing.Properties.VariableNames = ["Tempo","Any_Missing","Any_Consecutive","Async_Outliers","Tap_Miss","Tone_Miss","All_Miss","N_Events","NMiss_Tolerance","CritValue_Async","Conseq_Miss"];
    Results.ALL = [Results.File Results.BehavMeasures Results.AdaptParams Results.JointParams Results.Entropy Results.Missing];
    
    writetable(Results.ALL, 'C:/Users/marc1/OneDrive/Desktop/School/Aarhus University/Current Classes/Thesis/Analysis/ThesisADAM_boundedresults.csv');

    clearvars MissedCrit ConsecMiss CritAsync


    