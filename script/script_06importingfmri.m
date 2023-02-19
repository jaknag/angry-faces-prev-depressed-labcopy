%% Fighting with the dataset

clear;
close all;
% add fieldtrip and SPM paths
spm_path = '/Users/nagrodzkij/data/instal/spm12';
data_path = '/Users/nagrodzkij/data/';
output_path = '/Users/nagrodzkij/data/angry/output/';
input_path = '/Users/nagrodzkij/data/angry/input/';
motionparams_path = '/Users/nagrodzkij/data/angry/input/motion/';
scans_path = '/Users/nagrodzkij/data/angry/input/scans/';

cd(data_path);
if ~exist("fmri",'dir')
    mkdir fmri;
end

addpath(spm_path);

cd '/Users/nagrodzkij/data/angry/input/';
master = readtable('280_MRIIDs_03-06-15.csv');
load data_emoFace_excl.mat;
excluded_trials = readtable('excluded_trials.csv');

%% Make required folders in the fmri directory 
cd(output_path);
if ~exist('fmri','dir')
    mkdir('fmri');
end

fmri_path = sprintf('%sfmri/',output_path);

cd(fmri_path);
for i = 1:length(ccid_all)
    cd(fmri_path);
    ccid = ccid_all(i);
    fol_name = sprintf('sub-%.0f/',ccid);
    if ~exist(fol_name,'dir')
        mkdir(fol_name);
        cd(fol_name);
        mkdir func;
        mkdir anat;
        mkdir stats;
        mkdir sots;
    end
end

%% Prepare multiple condition files .mat FOR THREE CONDITIONS ONLY 
for iMultiple = 1:length(ccid_all_excl)
    ccid = ccid_all_excl(iMultiple);
    data = data_emoFace_excl(data_emoFace_excl.subj_idx == ccid,:);
    data_excluded = excluded_trials(excluded_trials.subj_idx == ccid,:);

    clearvars names onsets durations;

    names = cell(1,3);
    onsets = cell(1,3);
    durations = cell(1,3);
    
    ccid_str = string(ccid);
    particip_path = sprintf('%s/sub-%s/',fmri_path,ccid_str);
    sots_path = sprintf('%s/sots/',particip_path);

    for m = 1:2 % Loop through each of the four conditions
        conds = {'a_Angry','n_Neutral'};
        cond = conds{m};
        name = cond(3:end);
        emo = cond(1); % Get the emotion (a / n = Angry / Neutral)
        
        if strcmp(emo,'a') % Extract data for angry
            data_filtered = data(strcmp(data.condition,'Angry'),:);
           
        elseif strcmp(emo,'n') % Extract data for neutral
            data_filtered = data(strcmp(data.condition,'Neutral'),:);
        end

        onset = data_filtered.sot/1000; % Convert sots and button presses to seconds
        names{m} = name;
        onsets{m} = onset;
        durations{m} = 0;

    end

    names{3} = 'invalid';
    onsets{3} = data_excluded.sot/1000;
    if isempty(onsets{3})  % for participants who did not have any invalid trials
        onsets{3} = 0;
        fprintf('Participant %s did not have any invalid trials', ccid_str);
        sprintf('\n')
    end

    durations{3} = 0;

    condition_filename = sprintf('CC%s_condSpec_zero_3conds.mat',ccid_str);
    full_filename = [sots_path condition_filename];
        
    save(full_filename,"names","onsets","durations");
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Issues with filenames which needed to be resolved

% 1. Behavioural data 
% 
% Participant exists as 520079 in behavioural — their real CCID is 520097
% 
% 2. Scans 
% 
% Scans with following CCIDs have no equivalent in behavioural data:
% 
% Incorrect CCID -> CBUID -> Correct CCID 
% 
% CC222562 -> CBU140540 -> CC222652
% CC410113 -> CBU130953 -> CC410390
% CC410175 -> CBU131127 -> CC410179
% CC510350 -> CBU131056 -> CC510355
% CC610932 -> CBU140487 -> CC610392


% Put each scan in its participant folder and each motion param in its participant folder
clear iCopy;
clear ccid_mparam;
clear ccid_scan;

mparams = dir(sprintf('%s/*.txt', motionparams_path));
scans = dir(sprintf('%s/*.nii', scans_path));
cd '/Users/nagrodzkij/data/angry/output/';
save("scans_names","scans")

for iCopy = 1:numel(mparams)
    ccid_mparam(iCopy) = string(mparams(iCopy).name(3:8));
    motionparam_file = sprintf('%s/CC%s_motionParams.txt',motionparams_path,ccid_mparam(iCopy));
    particip_folder = sprintf('%s/sub-%s/func',fmri_path,ccid_mparam(iCopy));

    copyfile(motionparam_file,particip_folder);
end

for iCopy = 1:numel(scans)
    ccid_scan(iCopy) = string(scans(iCopy).name(17:22));
    scan_file = sprintf('%s/sswarfMR13010_CC%s-*.nii',scans_path,ccid_scan(iCopy));
    particip_folder = sprintf('%s/sub-%s/func',fmri_path,ccid_scan(iCopy));

    copyfile(scan_file,particip_folder);
end

