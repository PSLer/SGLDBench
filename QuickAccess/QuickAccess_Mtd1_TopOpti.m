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
MdlSelect = 'Bone'; %% Bone, Part, Part2, Part3, Bracket_GE, Molar, Fertility, Hanger, TopOptiShape
IO_LoadBuiltInDatasets(MdlSelect);
disp(['Preparing Voxel-based FEA Model Costs ', sprintf('%10.1f',toc(tStart)), 's'])

%%2. Optimization
DEBUG_ = 0; 
rMin_ = 2.5;
maxSharpness_ = 0.01;
nLoop_ = 30;
V_ = 0.4;
optimizer_ = 'OC';
constraintType_ = 'Global';
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(2, 0, 0);

TopOpti_CallTopOpti([]);

if ispc, system('"../src/quokka.exe" ../out/DesignVolume.nii'); end	
% if ispc, system('"../src/quokka.exe" ../out/alignmentMetricVolume_byStress.nii'); end