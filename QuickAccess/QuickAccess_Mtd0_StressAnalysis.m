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
MdlSelect = 'Bone'; %% Bone, Part, Part2, Part3, Bracket_GE, Molar, Fertility, Hanger, TopOptiShape
IO_LoadBuiltInDatasets(MdlSelect);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

% figure; view(gca,3);
% Vis_DrawMesh3D(gca, meshHierarchy_(1).boundaryNodeCoords, meshHierarchy_(1).boundaryEleFaces, 0);
% Vis_ShowLoadingCondition(gca, loadingCond_);
% Vis_ShowFixingCondition(gca, fixingCond_);

%% 2. Setup FEA
tStart = tic;
FEA_ApplyBoundaryCondition();
FEA_SetupVoxelBased();
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
