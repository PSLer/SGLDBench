clear all; clc;

%%Take care the routines when running on Linux
addpath('./src/');
addpath('./src/MEXfuncs/');
Data_GlobalVariables;
outPath_ = './out/';
inputVoxelfileName = './data/femur_B_R512.TopVoxel';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end
MEXfunc_ = true;
numLevels_ = []; %% [] Default
nonDyadic_ = 1; %%True or False
%%Data Loading
tStart = tic;
IO_ImportTopVoxels(inputVoxelfileName);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Settings
DEBUG_ = 0; 
constraintType_ = 'Local';
rMin_ = 2.6;
nLoop_ = 500;
maxSharpness_ = 0.01;
minChange_ = 1.0e-5;
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(2, 0, 0);
V_ = 0.53;
switch constraintType_
	case 'Global'		
		optimizer_ = 'OC';
	case 'Local'
		rHatMin_ = 8;
		alphaMin_ = 0.53;
end

%% Run
TopOpti_CallTopOpti([])
%% Less important
% profile off;
% profile viewer;
%%Vis.
% system('"./src/vape4d.exe" ./out/DesignVolume.nii');