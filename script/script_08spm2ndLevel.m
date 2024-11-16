%-----------------------------------------------------------------------
% Job saved on 03-Feb-2022 15:42:51 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

%% Set up paths 
clear;
close all;
% add fieldtrip and SPM paths
spm_path = '/Users/nagrodzkij/data/instal/spm12';
data_path = '/Users/nagrodzkij/data/';
output_path = '/Users/nagrodzkij/data/angry/output/';
input_path = '/Users/nagrodzkij/data/angry/input/';
motionparams_path = '/Users/nagrodzkij/data/angry/input/motion/';
scans_path = '/Users/nagrodzkij/data/angry/input/scans/';

model_dir = sprintf('%sfmri/2ndLevel/main/',output_path);

cd(output_path);
demog =readtable("table_demog_withedu_byccid.csv");

mask = sprintf('%ssdartel_GM_mask.nii.csv',input_path);
params_v = readtable(sprintf('%sparams_v_byparticip.csv',output_path));

%% Obtain the list of scans and contrasts 
cd '/Users/nagrodzkij/data/angry/output/';
scans = load("scans_names.mat");
scans = scans.scans;
all_contrasts = cell(numel(scans),1);

ccid_scan = nan(numel(scans),1);
for i = 1:numel(scans)
    ccid_scan(i) = string(scans(i).name(17:22));
    ccid_str = num2str(ccid_scan(i));
    participant_path = sprintf('%sfmri/sub-%s/',output_path,ccid_str);
    scan_filename = sprintf('%sstats/1stLevel/spmT_0004.nii',participant_path); % SPECIFY CONTRAST HERE
    all_contrasts{i} = (scan_filename);
end

IndexC = strfind(all_contrasts,'520055');
Index = find(not(cellfun('isempty',IndexC))); 
all_contrasts(Index,:)=[''];

ccid_scan = ccid_scan(ccid_scan~=520055);

all_contrasts_cellstr = cellstr(all_contrasts);

%% Filter demographics table by CCIDs for which we have scans 

[ccid_shared ccid_scan_idx ccid_demog_idx] = intersect(ccid_scan,demog.ccid);
demog_filtered = demog(ccid_demog_idx,:);

[ccid_shared ccid_scan_idx ccid_paramsv_idx] = intersect(ccid_scan,params_v.Var1);
params_v_filtered = params_v(ccid_paramsv_idx,:);

%% Prepare covariates 

% Load v_angry and v_neutral 
v_angry = params_v_filtered.v_angry;
v_neutral = params_v_filtered.v_neutral;

%%
matlabbatch{1}.spm.stats.factorial_design.dir = {model_dir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = all_contrasts_cellstr;

% Specify covariates 

matlabbatch{1}.spm.stats.factorial_design.cov(1).c = v_angry;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'v_angry';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1; % iCC = 1, centering around mean
% 
matlabbatch{1}.spm.stats.factorial_design.cov(2).c = v_neutral;
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'v_neutral';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1; 

matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Model estimation here
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrast manager 
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'faceF'; % F contrast for effect of age
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1];
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'face+'; % T contrast for positive effect of NEIS
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'face-'; % T contrast for negative effect of NEIS
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';


matlabbatch{3}.spm.stats.con.consess{4}.fcon.name = 'v_angryF'; % F contrast for effect of age
matlabbatch{3}.spm.stats.con.consess{4}.fcon.weights = [0 1];
matlabbatch{3}.spm.stats.con.consess{4}.fcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'v_angry+'; % T contrast for positive effect of NEIS
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 1];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'v_angry-'; % T contrast for negative effect of NEIS
matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 -1];
matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{7}.fcon.name = 'v_neutralF'; % F contrast for effect of age
matlabbatch{3}.spm.stats.con.consess{7}.fcon.weights = [0 0 1];
matlabbatch{3}.spm.stats.con.consess{7}.fcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'v_neutral+'; % T contrast for positive effect of NEIS
matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 1];
matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'v_neutral-'; % T contrast for negative effect of NEIS
matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [0 0 -1];
matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.delete = 0;


spm_jobman('run',matlabbatch);
