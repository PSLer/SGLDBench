clear all; clc;
addpath('./src/');
addpath('./src/MEXfuncs/');

%% Voxel Domain
tStart = tic;
Data_GlobalVariables;
inputVoxelfileName = './data/femur_B_R512.TopVoxel';
IO_ImportTopVoxels(inputVoxelfileName);
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(2, 0, 0);

%% Density Layout
fileDir = 'D:\wSpace\2024_pp_Summary3D\ressults\femur_R512_B\VF04\TopOpt_G\OC\';
numIts = 50;
for ii=0:1:numIts
	ii
	densityLayout_ = niftiread(strcat(fileDir, sprintf('intermeidateDensityLayout-It-%d.nii', ii)));
	densityLayout_ = flip(densityLayout_,1);
	densityLayout_ = densityLayout_(:);
	densityLayout_ = densityLayout_(meshHierarchy_(1).eleMapBack,1);
	densityLayout_(-1==densityLayout_) = 1;
	[~,minVal] = min(densityLayout_); densityLayout_(minVal) = 0;
	IO_ExportDesignInVolume_nii(strcat(fileDir, sprintf('new\intermeidateDensityLayout-It-%d_new.nii', ii)))
end

