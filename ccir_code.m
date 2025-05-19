clc;
subjstart=0;
subjend=1;
ch_offset=0;
nchan=7;
nchanr=1;  % number of channels per row. each channel has 2 graphs
alpha=[8,12];
theta=[4,8];
beta=[12,20];
delta=[0,4];
all=[0,20];
freq=all;
avgpow=false;
reassign=false;
idletime=[0,1];
% Initialize figure
if(~avgpow)
    figure;
end

% Define parameters
recording_duration = 60; % Total recording duration in seconds
fs = 500; % Sampling frequency in Hz (modify as needed)

% Define channel-to-brain region mapping
channel_regions = {
    'Fp1', 'Prefrontal Cortex';
    'Fp2', 'Prefrontal Cortex';
    'F3', 'Frontal Cortex';
    'F4', 'Frontal Cortex';
    'F7', 'Frontal Cortex';
    'F8', 'Frontal Cortex';
    'T3', 'Temporal Cortex';
    'T4', 'Temporal Cortex';
    'T5', 'Temporal Cortex';
    'T6', 'Temporal Cortex';
    'C3', 'Central Region';
    'C4', 'Central Region';
    'P3', 'Parietal Cortex';
    'P4', 'Parietal Cortex';
    'O1', 'Occipital Cortex';
    'O2', 'Occipital Cortex';
    'Fz', 'Midline (Frontal)';
    'Cz', 'Midline (Central)';
    'Pz', 'Midline (Parietal)';
    'A2_A1', 'Reference Electrode';
    'ECG', 'Heart Activity'
};
%creating tables after getting average power(look at end of code)
clear T;
clear T1;
clear T2;
varTypes = repelem({'string'}, 20);
varTypesPow = repelem({'double'}, 20);
varNames = channel_regions(1:20,1)';
sz=[36 20];
T  = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
T1 = table('Size',sz,'VariableTypes',varTypesPow,'VariableNames',varNames);
T2 = table('Size',sz,'VariableTypes',varTypesPow,'VariableNames',varNames);

% Path to folder containing EDF files
edf_folder = 'C:\Users\vidab\OneDrive\Desktop\eeg-mental-arithmetic\eeg-during-mental-arithmetic-tasks-1.0.0';

% Loop through subjects
% for subj = 0:35

for subj = subjstart:subjend %loop through subjects in a specific range
    % Define file names for EDF files
    subj_file_1 = fullfile(edf_folder, sprintf('Subject%02d_1.edf', subj));
    subj_file_2 = fullfile(edf_folder, sprintf('Subject%02d_2.edf', subj));

    % Check if files exist
    if exist(subj_file_1, 'file') && exist(subj_file_2, 'file')
        % Read EDF files
        subj_data_1 = edfread(subj_file_1);
        subj_data_2 = edfread(subj_file_2);
    else
        warning('EDF files for Subject%02d not found. Skipping.', subj);
        continue;
    end

    % Get the list of channel names from the field names
    channel_names = fieldnames(subj_data_1);

    % Filter out non-channel fields
    channel_names = channel_names(~ismember(channel_names, {'Properties', 'Record Time', 'Variables'}));
 
    %loop through channels
    for i = 1:nchan
       ch = i+ch_offset;
        % Extract channel name and region
        channel_name = channel_names{ch};
        % Remove EEG prefix
        channel_name_short = extractAfter(channel_name, 3);
        region_idx = find(strcmpi(channel_regions(:, 1), channel_name_short), 1);

        if isempty(region_idx)
            warning('Channel %s not found in brain region mapping.', channel_name);
            region = 'Unknown Region';
        else
            region = channel_regions{region_idx, 2};
        end

        % Extract data for the current channel
        x1 = subj_data_1.(channel_name);
        x2 = subj_data_2.(channel_name);
        % x1 is a 182 * 500 array. 500 samples mean 1 second worth of data
        % convert it into a linear array of 182*500=91000 entries
        newx1=[];
        for a=1:length(x1)
            newx1 = [newx1;x1{a}];
        end
        newx2=[];
        for a=1:length(x2)
            newx2 = [newx2;x2{a}];
        end
        t1 = (1:length(newx1)) / fs; % Generate time vector assuming uniform sampling
        t2 = (1:length(newx2)) / fs; % Generate time vector assuming uniform sampling

if(avgpow)
    %gathers average power
    [p, f] = pspectrum(newx1, fs, 'power');
    freqrange = [0 20]; 
    avgPower1 = bandpower(p, f, freqrange, 'psd');

else
    subplot(nchan/nchanr, nchanr*2, 2*i-1);
    %creates spectrogram
    pspectrum(newx1,fs,"spectrogram",'FrequencyLimits',freq,...
        OverlapPercent=0,MinThreshold=-20,TimeResolution=1,Reassign=reassign);
    caxis([-20,30]);
    xlim(idletime);
    title(['Before: ' channel_name_short])
end

if(avgpow)
    [p, f] = pspectrum(newx2, fs, 'power');
    freqrange = [0 20]; 
    avgPower2 = bandpower(p, f, freqrange, 'psd');

else
    subplot(nchan/nchanr, nchanr*2, 2*i);
    pspectrum(newx2,fs,"spectrogram",'FrequencyLimits',freq,...
        OverlapPercent=0,MinThreshold=-20,TimeResolution=1,Reassign=reassign);
    caxis([-20,30]);
    xlim([0,1]);
end

if(avgpow)
   if avgPower1>avgPower2
       final = 'lower'; %"lower" if average power of "before" is greater than average power of "during"
   end 
   if avgPower2>avgPower1 %"higher" if vice versa
       final = 'higher';
   end 
   if avgPower2 == avgPower1 %otherwise, its the same (unlikely)
       final = 'same';
   end
  T(subj+1,i) = {final};   
  T1(subj+1,i) = {avgPower1};
  T2(subj+1,i) = {avgPower2};
  %creates tables for results
end
    end
end
if(avgpow)
    writetable(T,'C:\Users\vidab\OneDrive\Documents\MATLAB\table.csv');
    writetable(T1,'C:\Users\vidab\OneDrive\Documents\MATLAB\table_before.csv');
    writetable(T2,'C:\Users\vidab\OneDrive\Documents\MATLAB\table_during.csv');
    %creates a csv file to import into Google Sheets
end
