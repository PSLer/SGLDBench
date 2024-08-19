function PSLs_Preparation4TSV()
	global outPath_;
	global meshHierarchy_;
	global boundingBox_;
	global nodeCoords_;
	global silhouetteStruct_;
	
	nx = meshHierarchy_(1).resX;
	ny = meshHierarchy_(1).resY;
	nz = meshHierarchy_(1).resZ;
	if 1
		nodeCoords_ = zeros(meshHierarchy_(1).numNodes,3);
		nodeCoords_(:,1) = double(niftiread(strcat(outPath_, 'cache_coordX.nii')));
		nodeCoords_(:,2) = double(niftiread(strcat(outPath_, 'cache_coordY.nii')));
		nodeCoords_(:,3) = double(niftiread(strcat(outPath_, 'cache_coordZ.nii')));
	else
		xSeed = boundingBox_(1,1):(boundingBox_(2,1)-boundingBox_(1,1))/nx:boundingBox_(2,1);
		ySeed = boundingBox_(2,2):(boundingBox_(1,2)-boundingBox_(2,2))/ny:boundingBox_(1,2);
		zSeed = boundingBox_(1,3):(boundingBox_(2,3)-boundingBox_(1,3))/nz:boundingBox_(2,3);
		nodeCoords_ = zeros(meshHierarchy_(1).numNodes,3);
		tmp = repmat(reshape(repmat(xSeed,ny+1,1), (nx+1)*(ny+1), 1), (nz+1), 1);
		nodPosX = reshape(tmp, ny+1, nx+1, nz+1);
		nodeCoords_(:,1) = tmp(meshHierarchy_(1).nodMapBack,1);
		tmp = repmat(repmat(ySeed,1,nx+1)', (nz+1), 1);
		nodPosY = reshape(tmp, ny+1, nx+1, nz+1);
		nodeCoords_(:,2) = tmp(meshHierarchy_(1).nodMapBack,1);
		tmp = reshape(repmat(zSeed,(nx+1)*(ny+1),1), (nx+1)*(ny+1)*(nz+1), 1);		
		nodPosZ = reshape(tmp, ny+1, nx+1, nz+1);
		nodeCoords_(:,3) = tmp(meshHierarchy_(1).nodMapBack,1);	
	end
	
	% iSurface = isosurface(nodPosX, nodPosY, nodPosZ, reshape(meshHierarchy_(1).nodMapForward, ny+1, nx+1, nz+1), 0);
	% iCap = isocaps(nodPosX, nodPosY, nodPosZ, reshape(meshHierarchy_(1).nodMapForward, ny+1, nx+1, nz+1), 0);
	% iCap.faces = size(iSurface.vertices,1) + iCap.faces;
	% silhouetteStruct_.vertices = [iSurface.vertices; iCap.vertices];
	% silhouetteStruct_.faces = [iSurface.faces; iCap.faces];
	% silhouetteStruct_ = CompactPatchVertices(silhouetteStruct_);	
end

function oPatchs = CompactPatchVertices(iPatchs)
	oPatchs = iPatchs;
	numOriginalVertices = size(iPatchs.vertices,1);
	numOriginalFaces = size(iPatchs.faces,1);
	validVertices = unique(iPatchs.faces);
	numValidVertices = size(validVertices,1);
	if numOriginalFaces==numOriginalVertices, return; end
	mapVerticesValid2Original = zeros(numOriginalVertices,1);
	mapVerticesValid2Original(validVertices) = (1:numValidVertices)';
	oPatchs.vertices = oPatchs.vertices(validVertices,:);
	oPatchs.faces = mapVerticesValid2Original(iPatchs.faces);
end