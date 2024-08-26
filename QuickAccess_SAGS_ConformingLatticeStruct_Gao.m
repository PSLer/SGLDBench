clear all; clc;

addpath('./src/');
addpath('./src/MEXfuncs/');
addpath('./external/Gao2017/');
addpath('./external/TetGen/');

%%Loading Voxel Model
tStart = tic;
Data_GlobalVariables;
inputVoxelfileName = './data/Voxel_R512.TopVoxel';
IO_ImportTopVoxels(inputVoxelfileName);
triSurfMeshfileName = './data/Tri_femur.ply';
IO_ImportSurfaceMesh(triSurfMeshfileName);
disp(['Prepare Voxel Model Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Setup for Simulation
tStart = tic;
if isempty(F_), FEA_ApplyBoundaryCondition(); end
if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
disp(['Setup FEA Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Stress Analysis on Solid Domain
tStart = tic;
densityLayout_ = ones(meshHierarchy_(1).numElements,1);
meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityLayout_(:));
Solving_AssembleFEAstencil();
U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
disp(['FEA Costs (Solid): ', sprintf('%10.3g',toc(tStart)) 's']);
ceList = TopOpti_ComputeUnitCompliance();
c0 = meshHierarchy_(1).eleModulus*ceList;
disp(['Compliance (Solid): ' sprintf('%10.5e ', c0)]);
tStart = tic;
dominantDirSolid = Common_ExtractDominantDirectionsFromPrincipalStressDirections();
disp(['Stress Analysis Costs (Solid): ', sprintf('%10.3g',toc(tStart)) 's']);

%% Stress Preparation for SAGS Modules
%% .1 Setup Getway Tet-mesh for Subsequent Input
SAGS_GenerateDelaunayTetMeshFromInputSurfaceMesh(10000);
%%.2 Interpolating Stress Field to Tet-mesh From High-res Cartesian Mesh
SAGS_InterpolatingStressFieldOnTetMesh();

%% Stress-aligned Conforming Lattice Infill Generation
edgeWidth = 2;
aspectRatio = 1.0;
targetDepositionRatio = 0.4;
numLayerboundary = 2;
numLayerLoads = 0;
numLayerFixation = 0;
SAGS_StressAlignedConformingLatticeGeneration(edgeWidth, targetDepositionRatio, numLayerboundary, ...
    numLayerLoads, numLayerFixation, aspectRatio);
%% Output Design
%%-1 for passive elements for vis
IO_ExportDesignInVolume_nii('./out/DesignVolume.nii');

%% Evaluating Compliance of Design
tStart = tic;
meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityLayout_);
Solving_AssembleFEAstencil();
maxIT_ = 500; tol_ = 0.01;
U_ = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
disp(['FEA Costs (Design): ', sprintf('%10.3g',toc(tStart)) 's']);
ceList = TopOpti_ComputeUnitCompliance();
c = meshHierarchy_(1).eleModulus*ceList;
disp(['Compliance (Design): ' sprintf('%10.5e ', c)]);

%% Stress Analysis on Solid Domain
tStart = tic;
dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections();
disp(['Stress Analysis (Cartesian Stress, von Mises Stress, Principal Stresses) Costs (Design): ', sprintf('%10.3g',toc(tStart)) 's']);

%% Compute Alignment Deviation
tStart = tic;
alignmentMetricVolume = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
niftiwrite(alignmentMetricVolume, './out/alignmentMetricVolume.nii');
disp(['Computing Stress Alignment Deviation Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% Visual Analytics
if 1
system('"./src/vape4d.exe" ./out/DesignVolume.nii');
else
system('"./src/vape4d.exe" ./out/alignmentMetricVolume.nii');
end