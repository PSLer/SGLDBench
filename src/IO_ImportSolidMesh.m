function IO_ImportSolidMesh(fileName)
	[~,~,dataType] = fileparts(fileName);
	switch dataType
		case '.mesh'
			IO_ImportSolidMesh_Format_mesh(fileName);
		case '.vtk'
			IO_ImportSolidMesh_Format_vtk(fileName);
		case '.msh'
			IO_ImportSolidMesh_Format_msh(fileName);
	end
end

function IO_ImportSolidMesh_Format_mesh(fileName)
	global inputSolidMesh_;
	global surfaceTriMesh_;
	
	%%Read Data
	fid = fopen(fileName, 'r');
	if ~strcmp(fscanf(fid, '%s', 1), 'MeshVersionFormatted'), error('Incompatible Mesh Data Format!'); end
	if 1~=fscanf(fid, '%d', 1), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'Dimension'), error('Incompatible Mesh Data Format!'); end
	if 3~=fscanf(fid, '%d', 1), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'Vertices'), error('Incompatible Mesh Data Format!'); end
	inputSolidMesh_.numNodes = fscanf(fid, '%d', 1);
	inputSolidMesh_.nodeCoords = fscanf(fid, '%f %f %f  %f', [4, inputSolidMesh_.numNodes])'; 
	inputSolidMesh_.nodeCoords(:,4) = [];
	meshType = fscanf(fid, '%s', 1);
	switch meshType
		case 'Hexahedra'
			inputSolidMesh_.meshType = 'HEX';
			inputSolidMesh_.numElements = fscanf(fid, '%d', 1);
			inputSolidMesh_.eNodMat = fscanf(fid, '%d %d %d %d %d %d %d %d %d', [9, inputSolidMesh_.numElements])'; 
			inputSolidMesh_.eNodMat(:,end) = [];
		case 'Tetrahedra'
			inputSolidMesh_.meshType = 'TET';
			inputSolidMesh_.numElements = fscanf(fid, '%d', 1)/4;
			inputSolidMesh_.eNodMat = fscanf(fid, '%d %d %d %d %d', [5, inputSolidMesh_.numElements])';
			inputSolidMesh_.eNodMat(:,end) = [];
		otherwise
			error('Incompatible Mesh Data Format!');	
	end
	fclose(fid);
	
	%%Extract Boundary Info of Mesh
	[inputSolidMesh_.boundaryNodeCoords, inputSolidMesh_.boundaryPatchNodMat, inputSolidMesh_.nodeState, ~] = ...
			Common_ExtractBoundaryInfoFromSolidMesh(inputSolidMesh_.nodeCoords, inputSolidMesh_.eNodMat);
	surfaceTriMesh_.nodeCoords = inputSolidMesh_.boundaryNodeCoords;
	switch inputSolidMesh_.meshType
		case 'HEX'
			surfaceTriMesh_.eNodMat = inputSolidMesh_.boundaryPatchNodMat(:,[1 2 3 3 4 1])';
			surfaceTriMesh_.eNodMat = reshape(surfaceTriMesh_.eNodMat(:), 3, numel(surfaceTriMesh_.eNodMat)/3)';
		case 'TET'
			surfaceTriMesh_.eNodMat = inputSolidMesh_.boundaryPatchNodMat;
	end
	surfaceTriMesh_.numNodes = size(surfaceTriMesh_.nodeCoords,1);
	surfaceTriMesh_.numElements = size(surfaceTriMesh_.eNodMat,1);
	inputSolidMesh_.state = 1;
	surfaceTriMesh_.state = 1;
end

function IO_ImportSolidMesh_Format_vtk(fileName)
	global inputSolidMesh_;
	global surfaceTriMesh_;
	
	%%Read Data
    fid = fopen(fileName, 'r');
	if ~strcmp(fscanf(fid, '%s', 1), '#'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'vtk'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'DataFile'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'Version'), error('Incompatible Mesh Data Format!'); end
	if 3.0~=fscanf(fid, '%f', 1), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'Volume'), error('Incompatible Mesh Data Format!'); end
    if ~strcmp(fscanf(fid, '%s', 1), 'mesh'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'ASCII'), error('Incompatible Mesh Data Format!'); end	
	if ~strcmp(fscanf(fid, '%s', 1), 'DATASET'), error('Incompatible Mesh Data Format!'); end	
	if ~strcmp(fscanf(fid, '%s', 1), 'UNSTRUCTURED_GRID'), error('Incompatible Mesh Data Format!'); end	
	if ~strcmp(fscanf(fid, '%s', 1), 'POINTS'), error('Incompatible Mesh Data Format!'); end
	inputSolidMesh_.numNodes = fscanf(fid, '%d', 1);
	if ~strcmp(fscanf(fid, '%s', 1), 'double'), error('Incompatible Mesh Data Format!'); end
	inputSolidMesh_.nodeCoords = fscanf(fid, '%f %f %f', [3, inputSolidMesh_.numNodes])';
	if ~strcmp(fscanf(fid, '%s', 1), 'CELLS'), error('Incompatible Mesh Data Format!'); end
	inputSolidMesh_.numElements = fscanf(fid, '%d', 1);
	if 9~=fscanf(fid, '%d', 1)/inputSolidMesh_.numElements, error('Incompatible Mesh Data Format!'); end
	inputSolidMesh_.eNodMat = fscanf(fid, '%d %d %d %d %d %d %d %d %d', [9, inputSolidMesh_.numElements])';
	inputSolidMesh_.eNodMat(:,1) = []; 
	inputSolidMesh_.eNodMat = inputSolidMesh_.eNodMat + 1;
	fclose(fid);
	
	%%Extract Boundary Info of Mesh
	[inputSolidMesh_.boundaryNodeCoords, inputSolidMesh_.boundaryPatchNodMat, inputSolidMesh_.nodeState, ~] = ...
			Common_ExtractBoundaryInfoFromSolidMesh(inputSolidMesh_.nodeCoords, inputSolidMesh_.eNodMat);
	surfaceTriMesh_.nodeCoords = inputSolidMesh_.boundaryNodeCoords;
	surfaceTriMesh_.eNodMat = inputSolidMesh_.boundaryPatchNodMat(:,[1 2 3 3 4 1])';
	surfaceTriMesh_.eNodMat = reshape(surfaceTriMesh_.eNodMat(:), 3, numel(surfaceTriMesh_.eNodMat)/3)';	
	inputSolidMesh_.meshType = 'HEX';
	surfaceTriMesh_.numNodes = size(surfaceTriMesh_.nodeCoords,1);
	surfaceTriMesh_.numElements = size(surfaceTriMesh_.eNodMat,1);
	inputSolidMesh_.state = 1;
	surfaceTriMesh_.state = 1;	
end

function IO_ImportSolidMesh_Format_msh(fileName)
	global inputSolidMesh_;
	global surfaceTriMesh_;
	
	%%Read Data
    fid = fopen(fileName, 'r');
	if ~strcmp(fscanf(fid, '%s', 1), '$MeshFormat'), error('Incompatible Mesh Data Format!'); end
	if 2.2~=fscanf(fid, '%f', 1), error('Incompatible Mesh Data Format!'); end
	if 0~=fscanf(fid, '%d', 1), error('Incompatible Mesh Data Format!'); end
	if 8~=fscanf(fid, '%d', 1), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), '$EndMeshFormat'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), '$Nodes'), error('Incompatible Mesh Data Format!'); end
	inputSolidMesh_.numNodes = fscanf(fid, '%d', 1);
	inputSolidMesh_.nodeCoords = fscanf(fid, '%d %f %f %f', [4, inputSolidMesh_.numNodes])';
	inputSolidMesh_.nodeCoords(:,1) = [];
	if ~strcmp(fscanf(fid, '%s', 1), '$EndNodes'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), '$Elements'), error('Incompatible Mesh Data Format!'); end
	inputSolidMesh_.numElements = fscanf(fid, '%d', 1);
	inputSolidMesh_.eNodMat = fscanf(fid, '%d %d %d %d %d %d %d', [7, inputSolidMesh_.numElements])';
	inputSolidMesh_.eNodMat(:,[1 2 3]) = [];
	fclose(fid);
	
	%%Extract Boundary Info of Mesh
	[inputSolidMesh_.boundaryNodeCoords, inputSolidMesh_.boundaryPatchNodMat, inputSolidMesh_.nodeState, ~] = ...
			Common_ExtractBoundaryInfoFromSolidMesh(inputSolidMesh_.nodeCoords, inputSolidMesh_.eNodMat);
	surfaceTriMesh_.nodeCoords = inputSolidMesh_.boundaryNodeCoords;
	surfaceTriMesh_.eNodMat = inputSolidMesh_.boundaryPatchNodMat;	
	inputSolidMesh_.meshType = 'TET';
	surfaceTriMesh_.numNodes = size(surfaceTriMesh_.nodeCoords,1);
	surfaceTriMesh_.numElements = size(surfaceTriMesh_.eNodMat,1);
	inputSolidMesh_.state = 1;
	surfaceTriMesh_.state = 1;		
end
