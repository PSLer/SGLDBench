clear all; clc;
addpath('./src/');
addpath('./src/MEXfuncs/');
addpath('./tempTest/');
% profile on;
%%Data Loading
tStart = tic;
Data_GlobalVariables;
inputVoxelfileName = './data/Voxel_R512.TopVoxel';
IO_ImportTopVoxels(inputVoxelfileName);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(3, 0, 0);

iPath = 'D:\wSpace\2024_pp_Summary3D\ressults\femur_R512\TopOpti_L\';
inputDensityLayoutFile = strcat(iPath, 'alignmentMetricVolume_byStress.nii');
outputDensityLayoutFile = strcat(iPath, 'alignmentMetricVolume_byStress_ExcludePassiveEles.nii');
densityLayout_ = niftiread(inputDensityLayoutFile);
densityLayout_ = flip(densityLayout_,1);
densityLayout_ = densityLayout_(:);
densityLayout_ = densityLayout_(meshHierarchy_(1).eleMapBack,1);
IO_ExportDesignInVolume_nii(outputDensityLayoutFile)

callRender = strcat('"./src/vape4d.exe"', char(strcat(" ", outputDensityLayoutFile)));
system(callRender);