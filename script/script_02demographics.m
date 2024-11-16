%% This script imports demographics, filters them by ccids which took part in the fMRI study
% And outputs a csv file with each experimental trial assigned various
% demographical data
%% Preparations 

clear;
close all;

spm_path = '/Users/nagrodzkij/data/instal/spm12';
fmri_path = '/Users/nagrodzkij/data/angry/input/scans/';
dataDir = '/Users/nagrodzkij/data/angry/input/beh/';
output_path = '/Users/nagrodzkij/data/angry/output/';
input_path = '/Users/nagrodzkij/data/angry/input/';

cd(output_path);
load('data_emoFace_excl.mat');
demog = load(sprintf('%s/demog/demog_genderDiscrimination.mat/',input_path));

%% Remove 'CC' from ccid 

if length(demog.demog.ccid{1})==8
    for k = 1 : length(demog.demog.ccid)
        cellContents = demog.demog.ccid{k};
        % Truncate and stick back into the cell
        demog.demog.ccid{k} = cellContents(3:end);
    end
else
    disp('CC has already been removed!');
end

%%
demog_table = struct2table(demog.demog);

ccid_all = ccid_all_excl;
ccid_demog = demog_table.ccid;
ccid_demog = str2double(ccid_demog);

%Find intersection between ccid_all_excl and the ccid in demographics, as
%well as their indices 
[common_ccid, ccid_all_idx, ccid_demog_idx] = intersect(ccid_all,ccid_demog);
%Filter the demographics table to obtain only the rows we need 
demog_filtered = demog_table(ccid_demog_idx,:);

%% Load accuracy and rt data and combine tables

accuracy_table = convertvars(accuracy_table,"ccid",'string');

demog_table = convertvars(demog_table,"ccid",'string');

table_accuracy_demog = join(accuracy_table,demog_table);
dob = table_accuracy_demog.dob;
dob = datetime(dob,"InputFormat","MM/uuuu");
age = years(session_date_excl-dob);
age_table = array2table(age);
demog_filtered = [demog_filtered,age_table];

table_accuracy_demog = removevars(table_accuracy_demog,"dob");
table_accuracy_demog = addvars(table_accuracy_demog, session_date_excl, dob, age,'After','ccid');

cd(output_path);
save("table_accuracy_demog","table_accuracy_demog");

%% Create a table where each trial is assigned age, hads_depression, hads_anxiety, acer. 

clearvars demog_bytrial;
for i = 1:height(data_emoFace_excl)
    idx_trial = string(data_emoFace_excl.subj_idx(i));
    demog_bytrial(i,:) = demog_filtered(strcmp(demog_filtered.ccid, idx_trial),:); %Filter the table by logical where index of a trial matches ccid
end 

%Expand the table to include accuracy, response time, stim col, response
%col
log_rt = log(data_emoFace_excl.rt);
accuracy_col = data_emoFace_excl.response == data_emoFace_excl.stim_col;

stim_gender(data_emoFace_excl.stim_col==0) = {'Male'};
stim_gender(data_emoFace_excl.stim_col==1) = {'Female'};
stim_gender = stim_gender';

demog_bytrial = addvars(demog_bytrial, stim_gender, ...
    data_emoFace_excl.response, accuracy_col, data_emoFace_excl.condition, data_emoFace_excl.rt, log_rt, ...
    'NewVariableNames',{'stim_gender','response','accuracy','emotion','rt','log_rt'});
%%
cd(output_path);
save("table_demog_bytrial","demog_bytrial");
save('demog_table','demog_filtered');

filename = 'table_demog_bytrial.csv';
writetable(demog_bytrial,filename);

filename = 'table_demog_byccid.csv';
writetable(demog_filtered,filename);