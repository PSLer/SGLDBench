clear all; clc;
if 0
addpath('./src/');
addpath('./src/MEXfuncs/');
addpath('./tempTest/');
% profile on;

%%Data Loading 

%% Voxel Domain
tStart = tic;
Data_GlobalVariables;
inputVoxelfileName = './data/Voxel_R512.TopVoxel';
IO_ImportTopVoxels(inputVoxelfileName);
%% Density Layout
densityLayout_ = niftiread('./out/Femur_R512/DesignVolume.nii'); densityLayout_ = flip(densityLayout_,1);
densityLayout_ = densityLayout_(:);
densityLayout_ = densityLayout_(meshHierarchy_(1).eleMapBack,1);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Setup
tStart = tic;
if isempty(F_), FEA_ApplyBoundaryCondition(); end
if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
disp(['Setup FEA Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Stress Analysis on Solid Domain
if 0
	tStart = tic;
	densityField = ones(meshHierarchy_(1).numElements,1);
	meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityField(:));
	Solving_AssembleFEAstencil();
	U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
	disp(['FEA Costs (Solid): ', sprintf('%10.3g',toc(tStart)) 's']);
	tStart = tic;
	dominantDirSolid = Common_ExtractDominantDirectionsFromPrincipalStressDirections();
    disp(['Dominant Direction Extraction (inc. Cartesian Stress, von Mises Stress, Principal Stresses) Costs (Solid): ', sprintf('%10.3g',toc(tStart)) 's']);
else
	dominantDirSolid = niftiread('./out/Femur_R512/dominantDirSolid.nii');
end


%% Stress Analysis on Design
if 0
	tStart = tic;
	meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityLayout_);
	Solving_AssembleFEAstencil();
	maxIT_ = 500; tol_ = 0.01;
	U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
	disp(['FEA Costs (Design): ', sprintf('%10.3g',toc(tStart)) 's']);
	tStart = tic;
	dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections();
    disp(['Stress Analysis (Cartesian Stress, von Mises Stress, Principal Stresses) Costs (Design): ', sprintf('%10.3g',toc(tStart)) 's']);
else
	dominantDirDesign = niftiread('./out/Femur_R512/dominantDirDesign.nii');
end

%% Compute Alignment Deviation
tStart = tic;
alignmentMetricVolume = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
niftiwrite(alignmentMetricVolume, './out/alignmentMetricVolume.nii');
disp(['Computing Stress Alignment Deviation Costs: ', sprintf('%10.3g',toc(tStart)) 's']);
end
%% Visual Analytics
system('"./src/vape4d.exe" D:\wSpace\2024_pp_Summary3D\ressults\femur_R1024\Voronoi\testR512\DesignVolume.nii');