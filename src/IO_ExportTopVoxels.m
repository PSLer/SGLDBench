function IO_ExportTopVoxels(fileName)
	global meshHierarchy_;
	global fixingCond_; 
	global loadingCond_;
	global passiveElements_;
	
	if 0==meshHierarchy_(1).state
		warning('Non Voxel Model is Available!'); return;
	end
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s %s %s %s', '#Voxel Model for SGLDBench'); fprintf(fid, '\n');
	fprintf(fid, '%s ', 'Resolution:'); 
	fprintf(fid, '%d %d %d\n', [meshHierarchy_(1).resX meshHierarchy_(1).resY meshHierarchy_(1).resZ]);
	fprintf(fid, '%s %s ', 'Solid voxels:');
	fprintf(fid, '%d\n', meshHierarchy_(1).numElements);
	fprintf(fid, '%d\n', meshHierarchy_(1).eleMapBack');
	fprintf(fid, '%s %s ', 'Passive elements:');
	fprintf(fid, '%d\n', numel(passiveElements_));
	if ~isempty(passiveElements_)
		fprintf(fid, '%d\n', meshHierarchy_(1).eleMapBack(passiveElements_));
	end	
	fprintf(fid, '%s ', 'Fixations:');
	fprintf(fid, '%d\n', size(fixingCond_,1));		
	if ~isempty(fixingCond_)
		fprintf(fid, '%d %d %d %d\n', [meshHierarchy_(1).nodMapBack(meshHierarchy_(1).nodesOnBoundary(fixingCond_(:,1))) fixingCond_(:,2:4)]');
	end
	fprintf(fid, '%s ', 'Loads:');
	fprintf(fid, '%d\n', size(loadingCond_,1));
	if ~isempty(loadingCond_)
		fprintf(fid, '%d %.4e %.4e %.4e\n', [double(meshHierarchy_(1).nodMapBack(meshHierarchy_(1).nodesOnBoundary(loadingCond_(:,1)))) loadingCond_(:,2:4)]');
	end
	fclose(fid);
end