clear all; clc;
addpath('./src/');
addpath('./src/MEXfuncs/');
addpath('./tempTest/');

%%Data Loading
tStart = tic;
Data_GlobalVariables;
inputVoxelfileName = './data/Voxel_R256.TopVoxel';
IO_ImportTopVoxels(inputVoxelfileName);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Settings
DEBUG_ = 0; 
constraintType_ = 'Global';
rMin_ = 1.6;
nLoop_ = 300;
maxSharpness_ = 0.1;
minChange_ = 1.0e-3;

switch constraintType_
	case 'Global'
		V_ = 0.4;
		optimizer_ = 'MMA';
	case 'Local'
		rHatMin_ = 18;
		alphaMin_ = 0.5;
end

%% Run
tStart = tic;
TopOpti_CallTopOpti([])
disp(['Run TopOpti Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
%% Less important

