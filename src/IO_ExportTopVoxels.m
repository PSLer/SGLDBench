function IO_ExportTopVoxels(fileName, varargin)
	opt = 2;
	switch opt
		case 1 %%Previous
			if 2==nargin
				IO_ExportTopVoxels_old(fileName, varargin{1});
			else
				IO_ExportTopVoxels_old(fileName);
			end		
		case 2 %%New
			if 2==nargin
				IO_ExportTopVoxels_new(fileName, varargin{1});
			else
				IO_ExportTopVoxels_new(fileName);
			end				
	end
end

function IO_ExportTopVoxels_new(fileName, varargin)
	global meshHierarchy_ fixingCond_ loadingCond_ passiveElements_ nelx_ nely_ nelz_;
	global densityLayout_;
	if 0==meshHierarchy_(1).state
		warning('Non Voxel Model is Available!'); return;
	end
	%%Squeeze Voxel Volume
	eleVolume = reshape(meshHierarchy_(1).eleMapForward, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	eleVolume(nely_+1:end,:,:) = [];
	eleVolume(:,nelx_+1:end,:) = [];
	eleVolume(:,:,nelz_+1:end) = [];
	solidVoxels = find(eleVolume);
	if ~isempty(passiveElements_)
		eleVolume = zeros(size(meshHierarchy_(1).eleMapForward));
		eleVolume(meshHierarchy_(1).eleMapBack(passiveElements_)) = 1;
		eleVolume = reshape(eleVolume, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
		eleVolume(nely_+1:end,:,:) = [];
		eleVolume(:,nelx_+1:end,:) = [];
		eleVolume(:,:,nelz_+1:end) = [];
		passiveEles = find(eleVolume);		
	end
	if ~isempty(fixingCond_)
		fixingCond_tmp = fixingCond_;
		fixingCond_tmp(:,1) = meshHierarchy_(1).nodesOnBoundary(fixingCond_(:,1));
		[~,sortAscending] = sort(fixingCond_tmp(:,1));
		fixingCond_tmp = fixingCond_tmp(sortAscending,:);
		nodVolume = zeros(size(meshHierarchy_(1).nodMapForward));
		nodVolume(meshHierarchy_(1).nodMapBack(fixingCond_tmp(:,1))) = 1;
		nodVolume = reshape(nodVolume, meshHierarchy_(1).resY+1, meshHierarchy_(1).resX+1, meshHierarchy_(1).resZ+1);
		nodVolume(nely_+2:end,:,:) = [];
		nodVolume(:,nelx_+2:end,:) = [];
		nodVolume(:,:,nelz_+2:end) = [];
		fixingCond_tmp(:,1) = find(nodVolume);
	end
	if ~isempty(loadingCond_)
		loadingCond_tmp = loadingCond_;
		loadingCond_tmp(:,1) = double(meshHierarchy_(1).nodesOnBoundary(loadingCond_tmp(:,1)));
		[~,sortAscending] = sort(loadingCond_tmp(:,1));
		loadingCond_tmp = loadingCond_tmp(sortAscending,:);
		nodVolume = zeros(size(meshHierarchy_(1).nodMapForward));
		nodVolume(meshHierarchy_(1).nodMapBack(loadingCond_tmp(:,1))) = 1;
		nodVolume = reshape(nodVolume, meshHierarchy_(1).resY+1, meshHierarchy_(1).resX+1, meshHierarchy_(1).resZ+1);
		nodVolume(nely_+2:end,:,:) = [];
		nodVolume(:,nelx_+2:end,:) = [];
		nodVolume(:,:,nelz_+2:end) = [];
		loadingCond_tmp(:,1) = double(find(nodVolume));		
	end	
	
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s %s %s %s', '#Voxel Model for SGLDBench'); fprintf(fid, '\n');
	fprintf(fid, '%s ', 'Version:'); fprintf(fid, '%.1f\n', 1.0);
	fprintf(fid, '%s ', 'Resolution:'); 	
	fprintf(fid, '%d %d %d\n', [nelx_ nely_ nelz_]);
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
		fprintf(fid, '%d %6.2e\n', [double(solidVoxels) densityLayout_]');
	else
		fprintf(fid, '%d\n', solidVoxels');
	end
	fprintf(fid, '%s %s ', 'Passive elements:');
	fprintf(fid, '%d\n', numel(passiveElements_));
	if ~isempty(passiveElements_)
		fprintf(fid, '%d\n', passiveEles);
	end	
	fprintf(fid, '%s ', 'Fixations:');
	fprintf(fid, '%d\n', size(fixingCond_,1));		
	if ~isempty(fixingCond_)
		fprintf(fid, '%d %d %d %d\n', fixingCond_tmp');
	end
	fprintf(fid, '%s ', 'Loads:');
	fprintf(fid, '%d\n', size(loadingCond_,1));
	if ~isempty(loadingCond_)
		fprintf(fid, '%d %.4e %.4e %.4e\n', loadingCond_tmp');
	end
	fprintf(fid, '%s %s ', 'Additional Loads:'); fprintf(fid, '%d\n', 0);	
	fclose(fid);
end

function IO_ExportTopVoxels_old(fileName, varargin)
	global meshHierarchy_ fixingCond_ loadingCond_ passiveElements_;
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