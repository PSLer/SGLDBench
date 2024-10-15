%% DEMO: PSLs-guided Infill Design
clear all; clc;
addpath('../');
addpath('../src/');
addpath('../src/MEXfuncs/');
Data_GlobalVariables;
outPath_ = '../out/';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end

%%1. Modeling
tStart = tic;
IO_ImportSurfaceMesh('../data/Tri_femur.ply');
MdlSelect = 'Bone'; %% Bone, Part, Part2, Part3, Bracket_GE, Molar, Fertility, Hanger, TopOptiShape
IO_LoadBuiltInDatasets(MdlSelect);
disp(['Preparing Voxel-based FEA Model Costs ', sprintf('%10.1f',toc(tStart)), 's'])

%%2. FEA && stress analysis
%%2.1 Building mesh hierarchy; Assembling computing stencil; Solving linear system
[complianceSolid_, ~] = FEA_ComputeComplianceVoxel();
%%2.2 Stress analysis
[cartesianStressField_, vonMisesStressField_] = FEA_StressAnalysis();  
dominantDirSolid = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressField_);

%%3. PSLs-guided lightweight design
%%3.1 Data Preparation
PSLs_Preparation4TSV();
%%3.2 Generation
majorPSLsOpt = true; mediumPSLsOpt = false; minorPSLsOpt = true; psDirIndicator = [majorPSLsOpt mediumPSLsOpt minorPSLsOpt];
V_ = 0.4; %%Prescribed material budget
edgeThickness = 3; %% #Layers of voxels around the PSL trajectories
passiveElesBoundary = 2; passiveElesLoads = 0; passiveElesFixations = 0;
PSLs_GeneratePSLsGuidedInfillDesign(psDirIndicator, edgeThickness, V_, passiveElesBoundary, passiveElesLoads, passiveElesFixations);
%%3.3 Output&Vis Design
fileName = strcat(outPath_, 'DesignVolume.nii');
IO_ExportDesignInVolume_Geo_nii(fileName);  
%%Show design with the local executable (Windows-only)
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
	%%Show alignment deviations with the local executable (Windows-only)
	% if ispc, system('"../src/quokka.exe" ../out/alignmentMetricVolume_byStress.nii'); end
	% if ispc, system('"../src/quokka.exe" ../out/alignmentMetricVolume_byEdge.nii'); end
end
