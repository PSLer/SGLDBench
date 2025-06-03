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
IO_LoadBuiltInDatasets('Bone');
disp(['Preparing Voxel-based FEA Model Costs ', sprintf('%10.1f',toc(tStart)), 's'])

%%2. Optimization
DEBUG_ = 0; 
rMin_ = 2.6;
maxSharpness_ = 0.01;
nLoop_ = 50;
V_ = 0.4;
optimizer_ = 'OC';
constraintType_ = 'Global';
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(2, 0, 0);

TopOpti_CallTopOpti([]);
% Open Web Renderer web('https://keksboter.github.io/quokka/') to investigate structural design './SGLDBench/out/ResultVolume_Design.nii'

% Optionally via the local executable renderer (Windows-only)
% if ispc, system('"../src/quokka_0-4-0.exe" ../out/ResultVolume_Design.nii'); end	
% if ispc, system('"../src/quokka_0-4-0.exe" ../out/ResultVolume_Design_StressAlignment.nii'); end
% if ispc, system('"../src/quokka_0-4-0.exe" ../out/ResultVolume_Design_vonMises.nii'); end
