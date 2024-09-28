clear all; clc;
addpath('../');
addpath('../src/');
addpath('../src/MEXfuncs/');

Data_GlobalVariables;
outPath_ = '../out/';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end

%%Data Loading 
%% Voxel Domain
tStart = tic;
inputVoxelfileName = '../data/Part_R512.TopVoxel';
IO_ImportTopVoxels(inputVoxelfileName);
%% Density Layout
densityLayout_ = niftiread('D:\wSpace\2024_pp_Summary3D\ressults\part_R512\TopOpti_G/DesignVolume.nii'); densityLayout_ = flip(densityLayout_,1);
densityLayout_ = densityLayout_(:);
densityLayout_ = densityLayout_(meshHierarchy_(1).eleMapBack,1);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Setup
tStart = tic;
if isempty(F_), FEA_ApplyBoundaryCondition(); end
if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
disp(['Setup FEA Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Stress Analysis on Solid Domain
tStart = tic;
densityField = ones(meshHierarchy_(1).numElements,1);
meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityField(:));
Solving_AssembleFEAstencil();
U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
disp(['FEA Costs (Solid): ', sprintf('%10.3g',toc(tStart)) 's']);
tStart = tic;
[cartesianStressField, ~] = FEA_StressAnalysis();
dominantDirSolid = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressField);
disp(['Dominant Direction Extraction (inc. Cartesian Stress, von Mises Stress, Principal Stresses) Costs (Solid): ', sprintf('%10.3g',toc(tStart)) 's']);

%% Stress Analysis on Design
tStart = tic;
meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityLayout_);
Solving_AssembleFEAstencil();
maxIT_ = 500; tol_ = 0.01;
U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
disp(['FEA Costs (Design): ', sprintf('%10.3g',toc(tStart)) 's']);
tStart = tic;
[cartesianStressField, ~] = FEA_StressAnalysis();
dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressField);
disp(['Stress Analysis (Cartesian Stress, von Mises Stress, Principal Stresses) Costs (Design): ', sprintf('%10.3g',toc(tStart)) 's']);

%% Compute Alignment Deviation
tStart = tic;
alignmentMetricVolume = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
niftiwrite(alignmentMetricVolume, '../out/alignmentMetricVolume_byStress.nii');
disp(['Computing Stress Alignment Deviation Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
%% Visual Analytics
% system('"./src/quokka.exe" ../out/alignmentMetricVolume_byStress.nii');