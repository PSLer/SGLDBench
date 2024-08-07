function FEA_SetupVoxelBased()
	global meshHierarchy_;
	global numLevels_;
	
	%%1. Initialize Element Stiffness Matrix
	meshHierarchy_(1).Ke = FEA_VoxelBasedElementStiffnessMatrix();	
	
	%%2. Initialize Solver
	%% Building Mesh Hierarchy for Geometric Multi-grid Solver
	Solving_BuildingMeshHierarchy();
		
	%% Transfer Boundary Condition to Coaser Levels
	fixedDOFsOnFiner = zeros(meshHierarchy_(1).numDOFs,1);
    fixedDOFsOnFiner(meshHierarchy_(1).fixedDOFs) = 1;
	for ii=2:numLevels_
		forcesOnCoarser = Solving_RestrictResidual(fixedDOFsOnFiner,ii);
if 0		
		meshHierarchy_(ii).freeDOFs = int32(find(0==forcesOnCoarser));
		meshHierarchy_(ii).fixedDOFs = setdiff((1:int32(meshHierarchy_(ii).numDOFs))', meshHierarchy_(ii).freeDOFs);
else
		meshHierarchy_(ii).freeDOFs = 0==forcesOnCoarser;
		meshHierarchy_(ii).fixedDOFs = true(size(meshHierarchy_(ii).freeDOFs));
		meshHierarchy_(ii).fixedDOFs(meshHierarchy_(ii).freeDOFs) = false;
end
		fixedDOFsOnFiner = forcesOnCoarser;
	end
	meshHierarchy_(1).Ks = meshHierarchy_(1).Ke;
end