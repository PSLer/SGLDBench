%% This script declares some global variables
% addpath('./external/Gao2017/');
%%Solid Mesh
global MEXfunc_; MEXfunc_ = true; %%MEX functionality is enabled to achieve faster indexing, one can disable it if there are any compiling issues, and Matlab's built-in indexing will be used
global inputSolidMesh_; inputSolidMesh_ = Data_ArbitraryMeshStruct();
global surfaceTriMesh_;	surfaceTriMesh_ = Data_ArbitraryMeshStruct();
global vertexEdgeGraph_; vertexEdgeGraph_ = Data_VertexEdgeGraphStruct();
global meshHierarchy_; meshHierarchy_ = Data_CartesianMeshStruct();
global boundingBox_; boundingBox_ = [];
%% Linear System Solver	
global tol_; tol_ = 1.0e-3; %% convergence tolerance of iterative linear system solver
global maxIT_; maxIT_ = 200; %% permitted maximum number of iteartion
global weightFactorJacobi_; weightFactorJacobi_ = 0.35; %% try reducing it in case failing to converge, scope: (0,1)
global coarsestResolutionControl_; coarsestResolutionControl_ = 10000;
global numLevels_;
global nonDyadic_; nonDyadic_ = 1; %%True or False
global optimizer_; optimizer_ = 'MMA'; %% 'MMA', 'OC'
global optimizerMovingStepSize_; optimizerMovingStepSize_ = 0.1;
%% Material properties
global modulus_; modulus_ = 1.0; %% Young's modulus	
global poissonRatio_; poissonRatio_ = 0.3;	%% Poisson's ratio
global modulusMin_; modulusMin_ = 1.0e-6;
global cellSize_; cellSize_ = 1;
global SIMPpenalty_; SIMPpenalty_ = 3;

%% Mesh & design domain description
% global surfTriMeshStruct_; %% input triangular surface mesh for voxelizing
global voxelizedVolume_; voxelizedVolume_ = []; %% voxelized model
% global characteristicSize_; characteristicSize_ = []; %% scalar or 'empty' (==the dimensionality of the bounding box of the input model)
global finestResolutionControl_; finestResolutionControl_ = 128; %% maximum number of elements along a single dimension
global nelx_; %% mesh resolution = nelx_ * nely_ * nelz_
global nely_;
global nelz_; 
% global boundingBox_;

%%GMGS parameters
%% minimum number of elements along a single dimension

global loadingCond_; loadingCond_ = []; %% applied forces
global fixingCond_; fixingCond_ = []; %% fixed nodes
global F_; F_ = []; %% force (right hand section)
global U_; U_ = []; %% displacement (solution of A*U_ = F_, A is system matrix)
global strainEnergyPerElement_; strainEnergyPerElement_ = [];
global cartesianStressField_; cartesianStressField_ = [];
global vonMisesStressField_; vonMisesStressField_ = [];
global complianceSolid_; complianceSolid_ = 0;
global complianceDesign_; complianceDesign_ = 0;
global volumeFractionDesign_; volumeFractionDesign_ = 0;
global voxelsOnBoundary_; voxelsOnBoundary_ = [];
global DEBUG_; DEBUG_ = 1;

% global iterationHist_; iterationHist_ = []; %% statistic of iterative solver
global densityLayout_; densityLayout_ = [];
global densityLayoutWithoutBoundary_; densityLayoutWithoutBoundary_ = [];
global passiveElements_; passiveElements_ = [];

global rMin_; rMin_ = 1.6;
global pMax_; pMax_ = 128;
global constraintType_; constraintType_ = 'Global';
global V_; V_ = 0.4;
global rHatMin_; rHatMin_ = 6;
global alphaMin_; alphaMin_ = 0.5;
global optimizer_; optimizer_ = 'MMA';
global move_; move_ = 0.1;
global nLoop_; nLoop_ = 300;
global maxSharpness_; maxSharpness_ = 0.05;
global minChange_; minChange_ = 1.0e-4;
global continueTopOpt_; continueTopOpt_ = 0;

%%penalty in SIMP
global penalty_; penalty_ = 3.0;
global penaltyIncrement_; penaltyIncrement_ = 0.2;
global penaltyUpdateIterations_; penaltyUpdateIterations_ = 30;
global p_; p_ = 16; 
global beta_; beta_ = 1;		% _beta, beta continuation
global eta_; eta_ = 0.5;		% projection threshold

%% node-pick operations (to facilitate applying for boundary condition)
global hdPickedNode_; hdPickedNode_ = [];
global hdSelectionBox_; hdSelectionBox_ = [];
global pickedNodeCache_; pickedNodeCache_ = [];

%%IO
global outPath_; outPath_ = './out/';



