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
if 1
	IO_ImportSurfaceMesh('../data/Tri_femur.ply');
	FEA_CreateVoxelizedModel(512);
	FEA_VoxelBasedDiscretization();
	loadingCond_ = load('../data/femur_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
	fixingCond_ = load('../data/femur_R512_fixa.bc');
else
    %%Use the IO functionalities in GUI to export *.TopVoxel File for repeated use
	IO_ImportTopVoxels('../data/Bearing_R512.TopVoxel'); %%Create from wrapped voxel file
end
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% 2. Setup FEA
tStart = tic;
if isempty(F_), FEA_ApplyBoundaryCondition(); end
if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
densityField = ones(meshHierarchy_(1).numElements,1); %%fully solid domain
meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityField(:));
disp(['Setup FEA Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

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
