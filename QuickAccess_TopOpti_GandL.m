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
disp(['Prepare Voxel Model Costs: ', sprinretf('%10.3g',toc(tStart)) 's']);

%% Settings
DEBUG_ = 0; 
constraintType_ = 'Local';
rMin_ = 1.6;
nLoop_ = 500;
maxSharpness_ = 0.01;
minChange_ = 1.0e-5;
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(2, 0, 0);
V_ = 0.6;
switch constraintType_
	case 'Global'		
		optimizer_ = 'OC';
	case 'Local'
		rHatMin_ = 12;
		alphaMin_ = 0.6;
end

%% Run
TopOpti_CallTopOpti([])
%% Less important
% profile off;
% profile viewer;
%%Vis.
system('"./src/vape4d.exe" ./out/DesignVolume.nii');