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
    %tStart = tic;
	global surfaceTriMesh_;

    fileContent = fileread(fileName);
    lines = strsplit(fileContent, '\n');

    vCount = sum(startsWith(lines, 'v '));
    fCount = sum(startsWith(lines, 'f '));
    nodeCoords = zeros(vCount, 3);
    eNodMat = zeros(fCount, 3);

    vIndex = 1;
    fIndex = 1;
    for i = 1:length(lines)
        line = lines{i};
        if startsWith(line, 'v ')
            nodeCoords(vIndex, :) = sscanf(line(3:end), '%f %f %f')';
            vIndex = vIndex + 1;
        elseif startsWith(line, 'f ')
            eNodMat(fIndex, :) = sscanf(line(3:end), '%d %d %d')';
            fIndex = fIndex + 1;
        end
    end

	surfaceTriMesh_.nodeCoords = nodeCoords;
	surfaceTriMesh_.eNodMat = eNodMat;
	surfaceTriMesh_.numNodes = size(surfaceTriMesh_.nodeCoords,1); 
	surfaceTriMesh_.numElements = size(surfaceTriMesh_.eNodMat,1);
	surfaceTriMesh_.state = 1;
    %fprintf("v: %d\n", surfaceTriMesh_.numNodes);
    %fprintf("f: %d\n", surfaceTriMesh_.numElements);
    %fprintf('t: %fs\n', toc(tStart));
end

function IO_ImportSurfaceMesh_Format_stl(fileName)
	global surfaceTriMesh_;
	FV = stlread(fileName);
	if 1
		surfaceTriMesh_.nodeCoords = FV.Points;
		surfaceTriMesh_.eNodMat = FV.ConnectivityList;
	else %%Try this if the above not working
		surfaceTriMesh_.nodeCoords = FV.vertices;
		surfaceTriMesh_.eNodMat = FV.faces;	
	end
	surfaceTriMesh_.numNodes = size(surfaceTriMesh_.nodeCoords,1); 
	surfaceTriMesh_.numElements = size(surfaceTriMesh_.eNodMat,1);
	surfaceTriMesh_.state = 1;	
end