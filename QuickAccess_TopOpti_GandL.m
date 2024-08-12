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
TopOpti_SetPassiveElements(2, 0, 0);
switch constraintType_
	case 'Global'
		V_ = 0.4;
		optimizer_ = 'OC';
	case 'Local'
		rHatMin_ = 18;
		alphaMin_ = 0.5;
end

%% Run
TopOpti_CallTopOpti([])
%% Less important

%%Vis.
targetVolumeFile = './out/DesignVolume.nii';
if ~exist(targetVolumeFile, 'file')
    if isempty(densityLayout_), warning('No design is available!'); return; end
    IO_ExportDesignInVolume_nii(targetVolumeFile);
end

if 1 %%Temporary for single format
    valueInput = niftiread('./out/DesignVolume.nii');
    valueOutput = double(valueInput);
    niftiwrite(valueOutput, './out/DesignVolume_double.nii');
    system('"./src/vape4d.exe" ./out/DesignVolume_double.nii');
else
    system('"./src/vape4d.exe" ./out/DesignVolume.nii');
end