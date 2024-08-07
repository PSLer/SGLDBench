clear all; clc;
addpath('./src/');
addpath('./src/MEXfuncs/');
addpath('./tempTest/');

%%Data Loading
tStart = tic;
Data_GlobalVariables;
inputVoxelfileName = './data/Voxel_R512.TopVoxel';
IO_ImportTopVoxels(inputVoxelfileName);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Setup
tStart = tic;
if isempty(F_), FEA_ApplyBoundaryCondition(); end
if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
densityField = ones(meshHierarchy_(1).numElements,1);
meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityField(:));
disp(['Setup FEA Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
% return
%% Assemble Computing Stencil
tStart = tic;
Solving_AssembleFEAstencil();
disp(['Assemble Computing Stencil Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
return
%% Iterative Solver
tStart = tic;
if 0
U_ = Solving_PreconditionedConjugateGradientSolver_previous(@Solving_KbyU_MatrixFree, @Solving_Vcycle_previous, F_, tol_, maxIT_, 'printP_ON');
else
U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
end
disp(['Linera System Solver Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Compute compliance
tStart = tic;
c = CmptComplianceTemp();
disp(['Compute Compliance Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
disp(['Compliance in total (weighted): ' sprintf('%10.5e ', c)]);

%% Less important

