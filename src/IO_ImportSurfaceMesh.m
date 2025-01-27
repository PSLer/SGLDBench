function IO_ImportSurfaceMesh(fileName)
	[~,~,dataType] = fileparts(fileName);
	switch dataType
		case '.ply'
			IO_ImportSurfaceMesh_Format_ply(fileName);
		case '.obj'
			IO_ImportSurfaceMesh_Format_obj(fileName);
		case '.stl'
			IO_ImportSurfaceMesh_Format_stl(fileName);
	end
end

function IO_ImportSurfaceMesh_Format_ply(fileName)
	global surfaceTriMesh_;
	%%Read Data
	fid = fopen(fileName, 'r');
	if ~strcmp(fscanf(fid, '%s', 1), 'ply'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'format'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'ascii'), error('Incompatible Mesh Data Format!'); end
	if 1.0~=fscanf(fid, '%f', 1), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'element'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'vertex'), error('Incompatible Mesh Data Format!'); end
	surfaceTriMesh_.numNodes = fscanf(fid, '%d', 1);
	if ~strcmp(fscanf(fid, '%s', 1), 'property'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'float'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'x'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'property'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'float'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'y'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'property'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'float'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'z'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'element'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'face'), error('Incompatible Mesh Data Format!'); end
	surfaceTriMesh_.numElements = fscanf(fid, '%d', 1);
	if ~strcmp(fscanf(fid, '%s', 1), 'property'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'list'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'uchar'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'uint'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'vertex_indices'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'end_header'), error('Incompatible Mesh Data Format!'); end
	surfaceTriMesh_.nodeCoords = fscanf(fid, '%f %f %f', [3, surfaceTriMesh_.numNodes])';
	surfaceTriMesh_.eNodMat = fscanf(fid, '%d  %d %d %d', [4, surfaceTriMesh_.numElements])';
	surfaceTriMesh_.eNodMat(:,1) = []; surfaceTriMesh_.eNodMat = surfaceTriMesh_.eNodMat + 1;
	fclose(fid);
	surfaceTriMesh_.state = 1;
end

function IO_ImportSurfaceMesh_Format_obj(fileName)
	global surfaceTriMesh_;
	%%Read Data
	nodeCoords = []; eNodMat = [];
	fid = fopen(fileName, 'r');
	while 1
		tline = fgetl(fid);
		if ~ischar(tline),   break,   end  % exit at end of file 
		ln = sscanf(tline,'%s',1); % line type 
		switch ln
			case 'v' % graph vertexs
				nodeCoords(end+1,1:3) = sscanf(tline(2:end), '%e')';
			case 'f'
				eNodMat(end+1,1:3) = sscanf(tline(2:end), '%d')';
		end
	end
	fclose(fid);
	surfaceTriMesh_.nodeCoords = nodeCoords;
	surfaceTriMesh_.eNodMat = eNodMat;
	surfaceTriMesh_.numNodes = size(surfaceTriMesh_.nodeCoords,1); 
	surfaceTriMesh_.numElements = size(surfaceTriMesh_.eNodMat,1);
	surfaceTriMesh_.state = 1;	
end

function IO_ImportSurfaceMesh_Format_stl(fileName)
	global surfaceTriMesh_;
	FV = stlread(fileName);
	surfaceTriMesh_.nodeCoords = FV.vertices;
	surfaceTriMesh_.eNodMat = FV.faces;
	surfaceTriMesh_.numNodes = size(surfaceTriMesh_.nodeCoords,1); 
	surfaceTriMesh_.numElements = size(surfaceTriMesh_.eNodMat,1);
	surfaceTriMesh_.state = 1;	
end