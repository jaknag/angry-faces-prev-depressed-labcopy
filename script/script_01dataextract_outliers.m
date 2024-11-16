%%This code extracts data from txt files

%% Preparations 

clear;
close all;

spm_path = '/Users/nagrodzkij/data/instal/spm12';
fmri_path = '/Users/nagrodzkij/data/angry/input/scans/';
dataDir = '/Users/nagrodzkij/data/angry/input/beh/';
output_path = '/Users/nagrodzkij/data/angry/output/';

%% Make a matrix of cbuIDs 

files = dir(sprintf('%s/*.txt', dataDir));
cbuID = nan(size(files,1), 1); %creates a matrix of not-numbers (cbuIDs) of size N x 1 where N is the number of files in the directory

for iFiles = 1:size(files,1)
    cbuID(iFiles) = str2double(files(iFiles).name(23:28));
end

%%
ccid_all = nan(length(cbuID),1);
proportion_responded = nan(length(cbuID), 1);

%session_date_all = nan(length(cbuID),1);

for iSub = 1:size(cbuID,1)
    
    fid = fopen(fullfile(files(iSub).folder, files(iSub).name), 'rt');
    
    bytesin = fread(fid)';
    
    data_cell = textscan(fid, '%s');
    
    % Heuristic to differentiate encoding from E-Prime 1 and 2
    if sum(bytesin==0) > 0
        textin = native2unicode(bytesin,'unicode');
    else
        textin = char(bytesin);
    end
    
    % Replace "colon+space" with a single, unused delimiter
    replaced = 0;
    delimiter = {'#' '%' '~' '|' '^'};
    for d = 1:length(delimiter)
        if isempty(strfind(textin,delimiter{d}))
            delimtext = strrep(textin,': ',delimiter{d});
            replaced = 1;
            break
        end
    end
    if replaced == 0
        error('Could not find a temporary delimiter to replace '': ''\nTry adding a new one to the line ''delimiter = {...};''\n')
    end

    allout = textscan(delimtext,'%s%s','delimiter',delimiter{d}); %Outputs a table with two columns: one column cotains what was before the delimiter, the other what was after
    c1 = allout{1};
    c2 = allout{2};
    
    fclose(fid);
    
    %Extract emotion, gender, response time, response accuracy 

    ind_emotion = find(strncmp('Emotion', c1, 7)); %Finds where the first 7 characters of strings in c1 match 'Emotion' (strncmp) and returns the row numbers for those
%     ind_emotion(1)=[]; ind_emotion(end)=[];
    emotion_cell = c2(ind_emotion);

    
    image_cell = c2(ind_emotion-8); %This is the cell containing the identity of the image used
    response_cell = c2(ind_emotion-6); %This is the cell containing the value of the response 
    correct_response_cell = c2(ind_emotion-7);
    response_accuracy = strncmp(response_cell, correct_response_cell, 2);
    
    ind_rt = ind_emotion - 5; %extracts the row number for response time
    RT_cell = c2(ind_rt); %extracts the response time value
    
    sot_cell = c2(ind_emotion-4);
    button_cell = str2double(sot_cell) + str2double(RT_cell);
    duration_cell = ones(length(ind_emotion),1);
    
    %create 3 arrays of response times, response (male/female) and stimulus (male/female)
    rt = zeros(length(RT_cell),1); 
    response = nan(length(RT_cell),1);
    stim_col = nan(length(RT_cell),1);
    
    for iRT = 1:length(RT_cell)
        rt(iRT) = str2double(RT_cell{iRT})/1000;
        if image_cell{iRT}(1) == 'm'
            stim_col(iRT) = 0; % male coded as zero
        elseif image_cell{iRT}(1) == 'f'
            stim_col(iRT) = 1; % female coded as one
        else 
            error('something went wrong with the image cell array')
        end
        
        % stimulus coding of response
        if strncmpi(response_cell(iRT),'26',2)
            response(iRT) = 1;
        elseif strncmpi(response_cell(iRT),'28',2)
            response(iRT) = 0;
        else
            response(iRT) = NaN;
        end
        
        if cbuID(iSub) == 131120 % this subject pressed different buttons! 
            if strncmpi(response_cell(iRT),'14',2)
                response(iRT) = 1;
            elseif strncmpi(response_cell(iRT),'22',2)
                response(iRT) = 0;
            else
                response(iRT) = NaN;
            end
            
        end
        
    end
    
    proportion_responded(iSub) = sum(~isnan(response))/length(response); %proportion responded 
    
    session_date = c2(strncmp('SessionDate', c1,11)); 
    session_date = session_date{1};
    session_date_all(iSub) = datetime(session_date, 'InputFormat', 'MM-dd-yyyy');

    ccid = c2(strncmp('Subject', c1, 7)); ccid = ccid{1};
    ccid_all(iSub) = str2double(ccid);
    cbu_id = c2(strncmp('CCDate', c1, 6)); cbu_id = str2double(cbu_id{1});
    if cbu_id ~= cbuID(iSub)
        error('something not right with the cbu id')
    end
    
    % preparing variables for table
    subj_idx = repmat(str2double(ccid), length(rt),1);
    condition = emotion_cell;
    sot = str2double(sot_cell);
    duration = duration_cell;
    button = button_cell;
    
    if iSub == 1
        data_emoFace_raw = table(subj_idx, sot, duration, button, response, rt, stim_col, condition);
    else
        data_emoFace_raw = [data_emoFace_raw; table(subj_idx, sot, duration, button, response, rt, stim_col, condition)];
    end
end  

session_date = transpose(session_date_all);

%% Number all trials by participant 

for i = 1:height(ccid_all)
    ccid = ccid_all(i); % Get the first ccid
    idx = find(data_emoFace_raw.subj_idx == ccid); % Filter the table by each ccid in turn
    
    for i = 1:height(idx)
        data_emoFace_raw(idx(i),"trial") = num2cell(i);
    end
end

%%
cd(output_path);
save("data_emoFace_raw.mat","data_emoFace_raw","cbuID","ccid_all","session_date");

% At the end of this section, data_emoFace is generated and saved.
% This is a table of: subject ID, whether they pressed male/female (male 0, female 1), what
% the response time was, what the stimulus was (male/female) (male 0,
% female 1) and what the condition was (Angry or Neutral). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exclude ultra-short reaction times and failed trials (where response was NaN)
% Create a table of excluded trials

data_emoFace_excl = data_emoFace_raw(data_emoFace_raw.rt > 0.25, :);

excluded_trials = data_emoFace_raw(find(data_emoFace_raw.rt<=0.25 | isnan(data_emoFace_raw.rt)),["subj_idx","trial","sot"]);

rows = height(data_emoFace_raw);
rows_good = height(data_emoFace_excl);
rows_bad = height(excluded_trials);

if rows_good + rows_bad == rows
    disp('Exclusion based on short reaction times complete');
else 
    disp('Something went wrong');
end


%% Find the % of stimuli that a participant responded to; the % of correct responses; and the ratio of correct responses / responses made
% It also finds the mean, max and min of those

% Create empty matrices for the values for later - this will ensure they
% will be the correct size

perc_responded = zeros(length(ccid_all),1);
perc_correct = zeros(length(ccid_all),1);
rate_correct = zeros(length(ccid_all),1);
failed_trials = zeros(length(ccid_all),1);

for iCor = 1:size(ccid_all, 1)
    
    ind_subj = find(ccid_all(iCor) == data_emoFace_raw.subj_idx);
    ind_subj_rtsexcl = find(ccid_all(iCor) == data_emoFace_excl.subj_idx);
    
    response_cell = data_emoFace_raw.response(ind_subj);
    response_cell_rtsexcl = data_emoFace_excl.response(ind_subj_rtsexcl);

    stim_cell_rtsexcl = data_emoFace_excl.stim_col(ind_subj_rtsexcl);

    response_accuracy = response_cell_rtsexcl == stim_cell_rtsexcl;

    failed_trials(iCor) = (length(response_cell) - length(response_cell_rtsexcl))/length(response_cell); %proportion of all trials which were ultra-short rts
    perc_responded(iCor) = sum(~isnan(response_cell_rtsexcl))/length(response_cell_rtsexcl); %percentage responded per participant
    perc_correct(iCor) = sum(response_accuracy)/length(response_cell_rtsexcl); %percentage correct per participant
    rate_correct(iCor) = perc_correct(iCor) / perc_responded(iCor);

    for iSanity = 1:size(ccid_all,1)
        if perc_responded(iCor) < perc_correct(iCor)
            disp('Something went wrong as number of responses is less than number of correct responses')
            break
        end
    end
end

disp ('Accuracy data - done');

accuracy_table = table(ccid_all,perc_responded,perc_correct,rate_correct,failed_trials);
accuracy_table = renamevars(accuracy_table,"ccid_all","ccid");

mean_perc_responded = mean(perc_responded);
min_perc_responded = min(perc_responded);
max_perc_responded = max(perc_responded);

mean_perc_correct = mean(perc_correct);
min_perc_correct = min(perc_correct);
max_perc_correct = max(perc_correct);

mean_rate_correct = mean(rate_correct);
min_rate_correct = min(rate_correct);
max_rate_correct = max(rate_correct);

%Automatic outlier identification (cut offs):

clear outlier_all
outlier_all = zeros(length(ccid_all),1);

for iOutl = 1:length(ccid_all)

    %INSERT CUT OFFS HERE

    cut_off_response = 0.75;
    cut_off_correct = 0.6;
    cut_off_rate_correct = 0.6;
    cut_off_rate_failed = 0.5;

    %if or(perc_responded(iOutl) < cut_off_response, perc_correct(iOutl) < cut_off_correct, rate_correct(iOutl) < cut_off_rate_correct)
    
    if perc_responded(iOutl) < cut_off_response || perc_correct(iOutl) < cut_off_correct || rate_correct(iOutl) < cut_off_rate_correct || failed_trials(iOutl) > cut_off_rate_failed
        fprintf(['An outlier seems to be ', num2str(ccid_all(iOutl)), '. \n', num2str(failed_trials(iOutl)) , ' of their trials were failed.\nThey responded to ', ...
        num2str(perc_responded(iOutl)*100), ' percent of all stimuli and they responded to ', ...
        num2str(perc_correct(iOutl)*100), ' percent of all stimuli correctly. \nThe rate of correct responses was ', ...
        num2str(rate_correct(iOutl)), ' of stimuli responded to.'])
        
        outlier_all(iOutl) = ccid_all(iOutl);
    end
end
outlier_all = outlier_all(outlier_all~=0);

disp(newline)
disp('Outliers analysed')
disp(newline)

%Exclude outliers based on participant number

[m,n]=size(data_emoFace_excl); %Get number of rows and number of columns of the table 

ccid_all_excl = ccid_all;
session_date_excl = session_date;

for i = 1:length(outlier_all)
    ccid_outlier = outlier_all(i); 

    idx_outlier_ccid = find(ccid_all==ccid_outlier);
    idx_outlier = find(data_emoFace_excl.subj_idx==ccid_outlier); %Find numbers of rows in the table corresponding to the outlier

    excluded_trials = [excluded_trials; data_emoFace_excl(idx_outlier,["subj_idx","trial","sot"])]; %Add those trials to the list of excluded trials

    data_emoFace_excl(idx_outlier,:) = []; %Delete those rows 
    accuracy_table(idx_outlier_ccid,:) = [];
    ccid_all_excl(idx_outlier_ccid,:) = [];
    session_date_excl(idx_outlier_ccid,:) = [];
    
end

excluded_trials = sortrows(excluded_trials); %Sort the excluded table by ccid 

%Sanity check -- check if number of excluded rows + number of rows after
%exclusion equals the number of rows overall 

rows = height(data_emoFace_raw);
rows_good = height(data_emoFace_excl);
rows_bad = height(excluded_trials);

if rows_good + rows_bad == rows
    disp('Exclusion based on participant complete');
else 
    disp('Something went wrong');
    return
end

%%
figure() 
histogram(data_emoFace_excl.rt)

%%
cd(output_path);
save("data_emoFace_excl.mat","excluded_trials","data_emoFace_excl","ccid_all_excl","session_date_excl","excluded_trials","accuracy_table","ccid_all");

%Saves csv file with all trials following outlier exclusion
writetable(data_emoFace_excl,"data_emoFace_excl.csv");

%Saves csv file with all excluded trials 
writetable(excluded_trials,"excluded_trials.csv");
