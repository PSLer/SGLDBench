%% DEMO: Stress-aware Voronooi Diagram Infill Design
clear all; clc;
addpath('../');
addpath('../src/');
addpath('../src/MEXfuncs/');
addpath('../externalModules/GradedVoronoiDiagram/');

Data_GlobalVariables;
outPath_ = '../out/';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end

%%1. Modeling
tStart = tic;
IO_ImportSurfaceMesh('../data/Tri_femur.ply');
FEA_CreateVoxelizedModel(512);
FEA_VoxelBasedDiscretization();
loadingCond_ = load('../data/femur_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
fixingCond_ = load('../data/femur_R512_fixa.bc');
disp(['Preparing Voxel-based FEA Model Costs ', sprintf('%10.1f',toc(tStart)), 's'])

%%2. FEA && stress analysis
%%2.1 Building mesh hierarchy; Assembling computing stencil; Solving linear system
[complianceSolid_, ~] = FEA_ComputeComplianceVoxel();
%%2.2 Stress analysis
[cartesianStressField_, vonMisesStressField_] = FEA_StressAnalysis();  
dominantDirSolid = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressField_);

%%3. Infill design
%%3.1 Data Preparation
numTetrahedraInGatewayMesh = 50000;
SAGS_GenerateDelaunayTetMeshFromInputSurfaceMesh(round(numTetrahedraInGatewayMesh));                      
SAGS_InterpolatingStressFieldOnTetMesh();
%%3.2 Generation
aspectRatio = 0.25;
V_ = 0.4; %%Prescribed material budget
edgeThickness = 3; %% #Layers of voxels around the PSL trajectories
passiveElesBoundary = 2; passiveElesLoads = 0; passiveElesFixations = 0;
SAGS_StressAwareGradedVoronoiDiagramGeneration(edgeThickness, V_, passiveElesBoundary, passiveElesLoads, passiveElesFixations, aspectRatio);
%%3.3 Output&Vis Design
fileName = strcat(outPath_, 'DesignVolume.nii');
IO_ExportDesignInVolume_Geo_nii(fileName);  
%%Show design with the local executable (Windows-only)
% figure; Vis_DrawGraph3D(gca, vertexEdgeGraph_.nodeCoords, vertexEdgeGraph_.eNodMat); light;
if ispc, system('"../src/quokka.exe" ../out/DesignVolume.nii'); end	

%%4. Design Evaluation
if 0
	maxIT_ = 500;
	tol_ = 1.0e-2; %%A slightly increased residual threshold for CG better balance efficiency and precision
	[complianceDesign_, volumeFraction_] = FEA_ComputeComplianceVoxel(densityLayout_);	
	% Solving_CG_GMGS('printP_ON'); %% Re-start CG without assembling computing stencil for a better converged solution
	[cartesianStressFieldDesign, ~] = FEA_StressAnalysis();  
	dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressFieldDesign);

    alignmentMetricVolumeByStressAlignment = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
    niftiwrite(alignmentMetricVolumeByStressAlignment, strcat(outPath_, 'alignmentMetricVolume_byStress.nii'));            
    alignmentMetricVolumeByEdgeAlignment = Common_ComputeEdgeAlignmentDeviation(dominantDirDesign);
    niftiwrite(alignmentMetricVolumeByEdgeAlignment, strcat(outPath_, 'alignmentMetricVolume_byEdge.nii'));	
	%%Show alignment deviations with the local executable (Windows-only)
	% if ispc, system('"../src/quokka.exe" ../out/alignmentMetricVolume_byStress.nii'); end
	% if ispc, system('"../src/quokka.exe" ../out/alignmentMetricVolume_byEdge.nii'); end
end
