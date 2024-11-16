%% Set up paths 
clear;
close all;
% add fieldtrip and SPM paths
spm_path = 'C:/Users/Joachim/Desktop/bunny/spm12/spm12'
data_path = 'F:/angry';
output_path = sprintf('%s/output/',data_path);
input_path = sprintf('%s/input/',data_path);
motionparams_path = sprintf('%s/motion/',input_path);

files2read = dir(sprintf('%sCC*',motionparams_path));

%%%%
subID = nan(length(files2read),1);
meanFD = nan(length(files2read),1);  % Mean FD for each participant
radius = 50;  % Convert rotations to mm

for iSub = 1:length(files2read)
    
    subID(iSub) = str2double(files2read(iSub).name(3:8));
    motionParams = dlmread(fullfile(files2read(iSub).folder, files2read(iSub).name));
    
    % Compute the frame-wise displacement (FD) for each time point
    FD_frames = nan(size(motionParams, 1) - 1, 1);  % Initialize FD for each frame

    for t = 2:size(motionParams, 1)
        % Compute the absolute differences in translation (X, Y, Z)
        trans_diff = abs(motionParams(t, 1:3) - motionParams(t-1, 1:3));
        
        % Compute the absolute differences in rotation (pitch, roll, yaw)
        rot_diff = abs(motionParams(t, 4:6) - motionParams(t-1, 4:6)) * radius;
        
        % Compute FD for this frame
        FD_frames(t-1) = sum(trans_diff) + sum(rot_diff);
    end
    
    % Compute the mean FD across all frames for this participant
    meanFD(iSub) = mean(FD_frames, 'omitnan');
    
end

% Display mean FD for each participant
disp(meanFD);
clf
histogram(meanFD)

movedidx = find(meanFD>0.5)

files = {files2read.name};
movedfiles = files(1,movedidx)

movedccids = nan(length(movedfiles),1)


for i = 1 : length(movedfiles)
    movedccids(i) = str2double(movedfiles{1,i}(3:8))
end

writematrix(movedccids,sprintf('%s/headmoved_ccids.csv',output_path));
