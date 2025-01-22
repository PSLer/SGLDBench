function IO_ExportTopVoxels(fileName, varargin)
	global meshHierarchy_;
	global fixingCond_; 
	global loadingCond_;
	global passiveElements_;
	global densityLayout_;
	if 0==meshHierarchy_(1).state
		warning('Non Voxel Model is Available!'); return;
	end
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s %s %s %s', '#Voxel Model for SGLDBench'); fprintf(fid, '\n');
	fprintf(fid, '%s ', 'Version:'); fprintf(fid, '%.1f\n', 1.0);
	fprintf(fid, '%s ', 'Resolution:'); 	
	fprintf(fid, '%d %d %d\n', [meshHierarchy_(1).resX meshHierarchy_(1).resY meshHierarchy_(1).resZ]);
	fprintf(fid, '%s %s', 'Density values:');
	densityValuesIncludedOpt = 0;
	if 2==nargin
		if varargin{1} 
			densityValuesIncludedOpt = 1;
		end
		if isempty(densityLayout_)
			warning('No Density Layout Available! Only Exports Voxel Model!'); 
			densityValuesIncludedOpt = 0;
		end		
	end		

	fprintf(fid, '%d\n', densityValuesIncludedOpt);
	fprintf(fid, '%s %s ', 'Solid voxels:');
	fprintf(fid, '%d\n', meshHierarchy_(1).numElements);
	if densityValuesIncludedOpt
		fprintf(fid, '%d %6.2e\n', [double(meshHierarchy_(1).eleMapBack) densityLayout_]');
	else
		fprintf(fid, '%d\n', meshHierarchy_(1).eleMapBack');
	end
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
	fprintf(fid, '%s %s ', 'Additional Loads:'); fprintf(fid, '%d\n', 0);	
	fclose(fid);
end