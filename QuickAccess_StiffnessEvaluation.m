clear all; clc;
%%Take care the routines when running on Linux
addpath('./src/');
addpath('./src/MEXfuncs/');
Data_GlobalVariables;
outPath_ = './out/';
inputVoxelfileName = './data/cantiR640_CMAME.TopVoxel';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end
MEXfunc_ = true;

%%Data Loading
tStart = tic;
IO_ImportTopVoxels(inputVoxelfileName);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Setup
tStart = tic;
if isempty(F_), FEA_ApplyBoundaryCondition(); end
if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
densityField = ones(meshHierarchy_(1).numElements,1);
meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityField(:));
disp(['Setup FEA Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Assemble Computing Stencil
tStart = tic;
Solving_AssembleFEAstencil();
disp(['Assemble Computing Stencil Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Iterative Solver
tStart = tic;
U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
disp(['Liner System Solver Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Compute compliance
tStart = tic;
ceList = TopOpti_ComputeUnitCompliance();
c = meshHierarchy_(1).eleModulus*ceList;
disp(['Compute Compliance Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
disp(['Compliance in total (weighted): ' sprintf('%10.5e ', c)]);

%% Less important

