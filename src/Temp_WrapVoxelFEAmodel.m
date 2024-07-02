function Temp_WrapVoxelFEAmodel(fileName)
	global meshHierarchy_;
	global fixingCond_; global loadingCond_;

	fid = fopen(fileName, 'w');	
	fprintf(fid, '%s %s %s', 'domain type: 3D'); fprintf(fid, '\n');
	fprintf(fid, '%s ', 'resolution:'); 
	fprintf(fid, '%d %d %d\n', [meshHierarchy_(1).resX meshHierarchy_(1).resY meshHierarchy_(1).resZ]);
	fprintf(fid, '%s %s ', 'valid elements:');
	fprintf(fid, '%d\n', meshHierarchy_(1).numElements);
	fprintf(fid, '%d\n', meshHierarchy_(1).eleMapBack');
	fprintf(fid, '%s %s ', 'fixed position:');
	fprintf(fid, '%d\n', size(fixingCond_,1));		
	if ~isempty(fixingCond_)
		fprintf(fid, '%d %d %d %d\n', [meshHierarchy_(1).nodMapBack(meshHierarchy_(1).nodesOnBoundary(fixingCond_(:,1))) fixingCond_(:,2:end)]');
	end
	fprintf(fid, '%s %s ', 'loading condition:');
	fprintf(fid, '%d\n', size(loadingCond_,1));
	if ~isempty(loadingCond_)
		fprintf(fid, '%d %.4e %.4e %.4e\n', [double(meshHierarchy_(1).nodMapBack(meshHierarchy_(1).nodesOnBoundary(loadingCond_(:,1)))) loadingCond_(:,2:end)]');
	end
	fprintf(fid, '%s %s %s', 'additional boundary conditions:');		
	fprintf(fid, ' %d\n', 0);
	fclose(fid);
end