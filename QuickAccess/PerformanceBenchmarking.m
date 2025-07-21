%% DEMO: This is to conduct performance benchmarking.
%% Consider a cuboid design domain discretized into 500x250x250 hexahedral elements, which correspond to about 100 million DOFs
%% With the prescribed boundary conditions and control parameters, run Topology Optimization on this domain
%% Comparing the processing time of this code on different hardware
clear all; clc;
addpath('../');
addpath('../src/');
addpath('../src/MEXfuncs/');

Data_GlobalVariables;
outPath_ = '../out/';
if ~exist(outPath_, 'dir'), mkdir(outPath_); end

%%1. Modeling
tStart = tic;
Shape_BuiltInCuboid(1, 0.5, 0.5);
FEA_CreateVoxelizedModel(500);
FEA_VoxelBasedDiscretization();
FEA_BuiltInBoundaryConditions4CuboidDesignDomain('Cuboid - Cantilever 3');
disp(['Preparing Voxel-based FEA Model Costs ', sprintf('%10.1f',toc(tStart)), 's']);

%%Show Problem Description
% axHandle_ = axes;
% Vis_DrawMesh3D(axHandle_, meshHierarchy_(1).boundaryNodeCoords, meshHierarchy_(1).boundaryEleFaces, 0);
% Vis_ShowLoadingCondition(axHandle_, loadingCond_);
% Vis_ShowFixingCondition(axHandle_, fixingCond_);
% view(axHandle_, 3);

%%2. Optimization
DEBUG_ = 0; 
rMin_ = 3;
maxSharpness_ = 0.01;
nLoop_ = 50;
V_ = 0.3;
optimizer_ = 'OC';
constraintType_ = 'Global';
[voxelsOnBoundary_, ~, ~] = TopOpti_SetPassiveElements(0, 0, 0);

TopOpti_CallTopOpti([]);
% Open Web Renderer web('https://keksboter.github.io/quokka/') to investigate structural design './SGLDBench/out/ResultVolume_Design.nii'
