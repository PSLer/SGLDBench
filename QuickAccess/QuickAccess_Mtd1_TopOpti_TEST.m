%% DEMO: Topology Optimization
clear all; clc;
addpath('../');
addpath('../src/');
addpath('../src/MEXfuncs/');

Data_GlobalVariables;
outPath_ = '../out/';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end
coarsestResolutionControl_ = 130;
tol_ = 1.0e-5;
nonDyadic_ = 0;
weightFactorJacobi_ = 0.6;

%%1. Modeling
tStart = tic;
IO_ImportTopVoxels('../data/test_canti_R64.TopVoxel'); %%Create from wrapped voxel file
disp(['Preparing Voxel-based FEA Model Costs ', sprintf('%10.1f',toc(tStart)), 's'])

%%2. Optimization
DEBUG_ = 0; 
rMin_ = 2.4;
maxSharpness_ = 0.01;
nLoop_ = 100;
V_ = 0.12;
optimizer_ = 'OC';
constraintType_ = 'Global_NEW';
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(0, 0, 0);

TopOpti_CallTopOpti([]);

% if ispc, system('"../src/quokka.exe" ../out/DesignVolume.nii'); end	
% if ispc, system('"../src/quokka.exe" ../out/alignmentMetricVolume_byStress.nii'); end