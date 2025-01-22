%% DEMO: Stress Analysis
clear all; clc;
addpath('../');
addpath('../src/');
addpath('../src/MEXfuncs/');

Data_GlobalVariables;
outPath_ = '../out/';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end

%%1. Data Loading
tStart = tic;
if 0
	MdlSelect = 'Part'; %% Bone, Part, Part2, Part3, Bracket_GE, Molar, Fertility, Hanger, TopOptiShape
	IO_LoadBuiltInDatasets(MdlSelect);
else
    %%Use the IO functionalities in GUI to export *.TopVoxel File for repeated use
	IO_ImportTopVoxels('../data/Part_R512.TopVoxel'); %%Create from wrapped voxel file
end
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% 2. Setup FEA
tStart = tic;cl
if isempty(F_), FEA_ApplyBoundaryCondition(); end
if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
densityLayout_ = niftiread('D:\wSpace\2024_pp_Summary3D\ressults\part_R512\TopOpti_L\new4Revision/DesignVolume.nii'); %densityLayout_ = flip(densityLayout_,1);
densityLayout_ = densityLayout_(:);
densityLayout_(-1==densityLayout_) = 1;
densityLayout_ = densityLayout_(meshHierarchy_(1).eleMapBack,1);
% densityLayout_ = ones(meshHierarchy_(1).numElements,1); %%fully solid domain
meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityLayout_(:));
disp(['Setup FEA Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
IO_ExportTopVoxels(strcat(outPath_, 'DesignVolume.TopVoxel'), 1);
return;
%% 3. Assemble Computing Stencil
tStart = tic;
Solving_AssembleFEAstencil();
disp(['Assemble Computing Stencil Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% 4. Solving FEA Linear System via Conjugate Gradien Method
tStart = tic;
U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
disp(['Liner System Solver Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% 5. Compute compliance
tStart = tic;
ceList = TopOpti_ComputeUnitCompliance();
c = meshHierarchy_(1).eleModulus*ceList;
disp(['Compliance in total (weighted): ' sprintf('%10.5e ', c)]);
disp(['Compute Compliance Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% 6. Compute Stress Field
tStart = tic;
[cartesianStressField_, vonMisesStressField_] = FEA_StressAnalysis();  
disp(['Compute Stress Field Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
