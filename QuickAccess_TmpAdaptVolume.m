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
fileDir = 'D:\wSpace\2024_pp_Summary3D\ressults\femur_R512_B\VF022\TopOpti_L\';
densityLayout_ = niftiread(strcat(fileDir, 'DesignVolume.nii')); 
densityLayout_ = flip(densityLayout_,1);
densityLayout_ = densityLayout_(:);
densityLayout_ = densityLayout_(meshHierarchy_(1).eleMapBack,1);
densityLayout_(-1==densityLayout_) = 1;
IO_ExportDesignInVolume_nii(strcat(fileDir, 'DesignVolume_new.nii'))
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
