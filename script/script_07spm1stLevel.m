%-----------------------------------------------------------------------
% Job saved on 31-Jan-2022 10:25:57 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
clear;
close all;
% add fieldtrip and SPM paths
spm_path = '/Users/nagrodzkij/data/instal/spm12';
data_path = '/Users/nagrodzkij/data/';
output_path = '/Users/nagrodzkij/data/angry/output/';
input_path = '/Users/nagrodzkij/data/angry/input/';
motionparams_path = '/Users/nagrodzkij/data/angry/input/motion/';
scans_path = '/Users/nagrodzkij/data/angry/input/scans/';

cd(output_path);
load data_emoFace_excl.mat;
scans = load("scans_names.mat");
scans = scans.scans;

user = getenv('USER'); % Will return the username for OSX operating systems

for i = 1:numel(scans)
    ccid_scan(i) = string(scans(i).name(17:22));
    subject = num2str(ccid_scan(i));

    %Prepare paths and filenames 
    preprocessed_filename = sprintf('/Users/%s/data/angry/output/fmri/sub-%s/func/sswarfMR13010_CC%s-0011.nii', user, subject, subject);
    subject_path = ['/Users/' user '/data/angry/output/fmri/sub-' subject '/'];
    model_destination_path = [subject_path 'stats/1stLevel'];
    modelSpec_filename = [subject_path 'sots/CC' subject '_condSpec_zero_3conds.mat'];
    motionParams_filename = [subject_path '/func/CC' subject '_motionParams.txt'];

    %Check if relevant scan is available 
    if exist([preprocessed_filename]) == 0
        display('Scan not available')
    end    

    matlabbatch{1}.spm.stats.fmri_spec.dir = {model_destination_path};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 32;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 16;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {preprocessed_filename};
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {modelSpec_filename};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {motionParams_filename};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Average effect of valid faces';
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1 1];
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.consess{2}.fcon.name = 'Main effect of Emotion';
    matlabbatch{3}.spm.stats.con.consess{2}.fcon.weights = [1 -1];
    matlabbatch{3}.spm.stats.con.consess{2}.fcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.consess{3}.fcon.name = 'Main effect of validity';
    matlabbatch{3}.spm.stats.con.consess{3}.fcon.weights = [0.5 0.5 -1];
    matlabbatch{3}.spm.stats.con.consess{3}.fcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Positive effect of Emotion (Angry)';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1 -1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'Positive effect of validity';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0.5 0.5 -1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Positive effect of Angry vs null';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1 0];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Negative effect of Angry vs null';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [-1 0];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.consess{8}.fcon.name = 'Effects_of_interest';
    matlabbatch{3}.spm.stats.con.consess{8}.fcon.weights = [eye(2)];
    matlabbatch{3}.spm.stats.con.consess{8}.fcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.delete = 0;

    spm_jobman('run',matlabbatch);

end

