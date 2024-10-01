%% DEMO: Porous Infill Optimization
clear all; clc;
addpath('../');
addpath('../src/');
addpath('../src/MEXfuncs/');

Data_GlobalVariables;
outPath_ = '../out/';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end

%%Data Loading
tStart = tic;
if 1
	IO_ImportSurfaceMesh('../data/Tri_femur.ply');
	FEA_CreateVoxelizedModel(512);
	FEA_VoxelBasedDiscretization();
	loadingCond_ = load('../data/femur_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
	fixingCond_ = load('../data/femur_R512_fixa.bc');
end
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%%2. Optimization
DEBUG_ = 0; 
maxSharpness_ = 0.01;
nLoop_ = 500;
rMin_ = 2.6;
rHatMin_ = 8;
alphaMin_ = 0.53;
constraintType_ = 'Local';
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(2, 0, 0);

TopOpti_CallTopOpti([]);

if ispc, system('"../src/quokka.exe" ../out/DesignVolume.nii'); end	
% if ispc, system('"../src/quokka.exe" ../out/alignmentMetricVolume_byStress.nii'); end
