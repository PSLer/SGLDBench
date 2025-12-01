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
if 0
    MdlSelect = 'Bone'; %% Bone, Part, Part2, Part3, Bracket_GE, Molar, Fertility, Hanger, TopOptiShape
    IO_LoadBuiltInDatasets(MdlSelect);
else
    Shape_BuiltInCuboid(1, 0.5, 0.5);
    FEA_CreateVoxelizedModel(500);
    FEA_VoxelBasedDiscretization();
    FEA_BuiltInBoundaryConditions4CuboidDesignDomain('Cuboid - Cantilever 3');
end
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
[U_, tSpMV, tPtV, its] = Solving_PreconditionedConjugateGradientSolver_TMP(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, 'printP_ON');
tTotal = toc(tStart);
disp(['MGCG in total Costs: ', sprintf('%10.3g',tTotal) 's']);
disp(['MGCG per iteration Costs: ', sprintf('%10.3g',tTotal/its) 's']);
disp(['SpMV per iteration costs: ', sprintf('%10.3g',tSpMV/its) 's']);
disp(['PtV per iteration costs: ', sprintf('%10.3g',tPtV/its) 's']);
disp(['Others per iteration costs: ', sprintf('%10.3g',tTotal/its-tSpMV/its-tPtV/its) 's']);
%% 5. Compute compliance
% tStart = tic;
% ceList = TopOpti_ComputeUnitCompliance();
% c = meshHierarchy_(1).eleModulus*ceList;
% disp(['Compliance in total: ' sprintf('%10.5e ', c)]);
% disp(['Compute Compliance Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

%% 6. Compute Stress Field
% tStart = tic;
% [cartesianStressField_, vonMisesStressField_] = FEA_StressAnalysis();  
% disp(['Compute Stress Field Costs: ', sprintf('%10.3g',toc(tStart)) 's']);

function [y, tMtV_, tPtV_, varargout] = Solving_PreconditionedConjugateGradientSolver_TMP(AtX, PtV, b, tol, maxIT, printP, varargin)
	%%0. arguments introduction
	%%AtX --- function handle for the product of system matrix and vector
	%%b --- right hand section
	%%tol --- stopping condition: resnrm < discrepancy
	%%maxIT --- mAtXximum number of iterations
tMtV_ = 0;
tPtV_ = 0;
	normB = norm(b);
	if 7==nargin, y = varargin{1}; else, y = zeros(size(b)); end
tStart1 = tic;		
	rVec1 = b - AtX(y);
tMtV_ = tMtV_ + toc(tStart1);
tStart2 = tic;		
	zVec = PtV(rVec1);
tPtV_ = tPtV_ + toc(tStart2);		
	pVec = zVec;
	x1Val = zVec' * rVec1;
	for its=1:maxIT
tStart1 = tic;		
		zVec = AtX(pVec);
tMtV_ = tMtV_ + toc(tStart1);			
		lambda = x1Val / (pVec' * zVec);
		y = y + lambda * pVec;
		rVec1 = rVec1 - lambda*zVec;
		resnorm = norm(rVec1)/normB;
		if printP(1)
			disp([' It.: ' sprintf('%4i',its) ' Res.: ' sprintf('%16.6e',resnorm)]);
		end		
		if resnorm<tol
			if printP(2)
				disp(['CG solver converged at iteration' sprintf('%5i', its) ' to a solution with relative residual' sprintf('%16.6e',resnorm)]);
			end
			break;
		end
tStart2 = tic;			
		zVec = PtV(rVec1);
tPtV_ = tPtV_ + toc(tStart2);		
		x2Val = zVec' * rVec1;
		pVec = zVec + x2Val / x1Val * pVec;
		x1Val = x2Val;
	end

	if its == maxIT
		warning('Exceed the maximum iterate numbers');
		disp(['The iterative process stops at residual = ' sprintf('%10.4f',resnorm)]);		
	end
	if 4==nargout, varargout{1} = its; end
end
