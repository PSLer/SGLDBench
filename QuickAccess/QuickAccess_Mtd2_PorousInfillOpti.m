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
if 1 %%To save size of the code repository
	MdlSelect = 'Bone'; %% Bone, Part, Part2, Part3, Bracket_GE, Molar, Fertility, Hanger, TopOptiShape
	IO_LoadBuiltInDatasets(MdlSelect);
else %% This is for general use. 	
	%One can create the voxel model via the GUI, and export it in the tailored .TopVoxel format, then run the simulation with this script
	IO_ImportTopVoxels('../data/NAME.TopVoxel'); %%Create from wrapped voxel file
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

% if ispc, system('"../src/quokka_0-4-0.exe" ../out/ResultVolume_Design.nii'); end	
% if ispc, system('"../src/quokka_0-4-0.exe" ../out/ResultVolume_Design_StressAlignment.nii'); end
% if ispc, system('"../src/quokka_0-4-0.exe" ../out/ResultVolume_Design_vonMises.nii'); end
