function IO_ExportCoarserVoxelMdl(iLevel)
	global meshHierarchy_ F_;

	if 0==meshHierarchy_(1).state
		warning('Non Voxel Model is Available!'); return;
	end
	if iLevel<=1, return; end
	
	fixedDOFsOnFiner = zeros(meshHierarchy_(1).numDOFs,1);
    fixedDOFsOnFiner(meshHierarchy_(1).fixedDOFs) = 1;	
	forcesOnFiner = full(F_);
	
	for ii=2:iLevel
		fixedDOFsOnFiner = Solving_RestrictResidual(fixedDOFsOnFiner,ii);
		forcesOnFiner = Solving_RestrictResidual(forcesOnFiner,ii);
	end
	fixedDOFsOnFiner(fixedDOFsOnFiner>0) = 1;
	fixedDOFsOnFiner = reshape(fixedDOFsOnFiner, 3, meshHierarchy_(iLevel).numNodes)';
	forcesOnFiner = reshape(forcesOnFiner, 3, meshHierarchy_(iLevel).numNodes)';
	tmp_fixedDOFsOnFiner = abs(fixedDOFsOnFiner); tmp_fixedDOFsOnFiner = sum(tmp_fixedDOFsOnFiner,2);
	tmp_forcesOnFiner = abs(forcesOnFiner); tmp_forcesOnFiner = sum(tmp_forcesOnFiner,2);
	fixedNodesRaw = find(tmp_fixedDOFsOnFiner>0);
	forceNodesRaw = find(tmp_forcesOnFiner>0);
	fixedNodes = intersect(fixedNodesRaw, meshHierarchy_(iLevel).nodesOnBoundary);
	forceNodes = intersect(forceNodesRaw, meshHierarchy_(iLevel).nodesOnBoundary);
	fixingCond = [meshHierarchy_(iLevel).nodMapBack(fixedNodes) fixedDOFsOnFiner(fixedNodes,:)];
	loadingCond = [double(meshHierarchy_(iLevel).nodMapBack(forceNodes)) forcesOnFiner(forceNodes,:)];
	
	resX = meshHierarchy_(iLevel).resX;
	resY = meshHierarchy_(iLevel).resY;
	resZ = meshHierarchy_(iLevel).resZ;
	numElements = meshHierarchy_(iLevel).numElements;
	eleMapBack = meshHierarchy_(iLevel).eleMapBack;
	
	fileName = '../out/FEMdlcoarser.TopVoxel';
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s %s %s %s', '#Voxel Model for SGLDBench'); fprintf(fid, '\n');
	fprintf(fid, '%s ', 'Version:'); fprintf(fid, '%.1f\n', 1.0);
	fprintf(fid, '%s ', 'Resolution:'); 	
	fprintf(fid, '%d %d %d\n', [resX resY resZ]);
	fprintf(fid, '%s %s', 'Density values:');
	densityValuesIncludedOpt = 0;
	fprintf(fid, '%d\n', densityValuesIncludedOpt);
	fprintf(fid, '%s %s ', 'Solid voxels:');
	fprintf(fid, '%d\n', numElements);
	fprintf(fid, '%d\n', eleMapBack');
	fprintf(fid, '%s %s ', 'Passive elements:');
	fprintf(fid, '%d\n', 0);
	fprintf(fid, '%s ', 'Fixations:');
	fprintf(fid, '%d\n', size(fixingCond,1));		
	if ~isempty(fixingCond)
		fprintf(fid, '%d %d %d %d\n', fixingCond');
	end
	fprintf(fid, '%s ', 'Loads:');
	fprintf(fid, '%d\n', size(loadingCond,1));
	if ~isempty(loadingCond)
		fprintf(fid, '%d %.4e %.4e %.4e\n', loadingCond');
	end
	fprintf(fid, '%s %s ', 'Additional Loads:'); fprintf(fid, '%d\n', 0);	
	fclose(fid);
end