%% DEMO: Topology Optimization
clear all; clc;
addpath('../');
addpath('../src/');
addpath('../src/MEXfuncs/');

Data_GlobalVariables;
outPath_ = '../out/';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end

%%1. Modeling
tStart = tic;
if 0
	IO_ImportSurfaceMesh('../data/Tri_femur.ply');
	FEA_CreateVoxelizedModel(512);
	FEA_VoxelBasedDiscretization();
	loadingCond_ = load('../data/femur_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
	fixingCond_ = load('../data/femur_R512_fixa.bc');
else
	IO_ImportTopVoxels('../data/part_R256.TopVoxel'); %%Create from wrapped voxel file
end
disp(['Preparing Voxel-based FEA Model Costs ', sprintf('%10.1f',toc(tStart)), 's'])

%%2. Optimization
DEBUG_ = 0; 
rMin_ = 2.5;
maxSharpness_ = 0.01;
nLoop_ = 20;
V_ = 0.3;
optimizer_ = 'OC';
constraintType_ = 'Global';
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(0, 3, 3);

TopOpti_CallTopOpti([]);

if ispc, system('"../src/quokka.exe" ../out/DesignVolume.nii'); end	
% if ispc, system('"../src/quokka.exe" ../out/alignmentMetricVolume_byStress.nii'); end