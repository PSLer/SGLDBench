function Temp_IO_ExportCartesianGrid_Mesh()
	global outPath_;
	global meshHierarchy_;

	nodeCoords_ = zeros(meshHierarchy_(1).numNodes,3);
	nodeCoords_(:,1) = double(niftiread(strcat(outPath_, 'cache_coordX.nii')));
	nodeCoords_(:,2) = double(niftiread(strcat(outPath_, 'cache_coordY.nii')));
	nodeCoords_(:,3) = double(niftiread(strcat(outPath_, 'cache_coordZ.nii')));
	
	fileName = strcat(outPath_, 'CartesianGridHex.mesh');
	fid = fopen(fileName, 'w');	
	%%header
	fprintf(fid, '%s ', 'MeshVersionFormatted');
	fprintf(fid, '%d\n', 1);
	fprintf(fid, '%s ', 'Dimension'); 
	fprintf(fid, '%d\n', 3);
	fprintf(fid, '%s ', 'Vertices');
	fprintf(fid, '%d\n', meshHierarchy_(1).numNodes);
	%%nodes
	fprintf(fid, '%.6e %.6e %.6e %d\n', [nodeCoords_ zeros(meshHierarchy_(1).numNodes,1)]');
	%%Cells
	fprintf(fid, '%s ', 'Hexahedra');
	fprintf(fid, '%d \n', meshHierarchy_(1).numElements);
	fprintf(fid, '%d %d %d %d %d %d %d %d %d\n', [meshHierarchy_(1).eNodMat zeros(meshHierarchy_(1).numElements,1)]');
	fprintf(fid, '%s', 'End');
	fclose(fid);

	edgeGraph = Data_VertexEdgeGraphStruct();
	cell2edge = [1 2  2 3  3 4  4 1  5 6  6 7  7 8  8 5  1 5  2 6  3 7  4 8];
	numEdgePerCell = 12;
	edgeIndices = meshHierarchy_(1).eNodMat(:, cell2edge)';
	edgeIndices = reshape(edgeIndices, 2, numEdgePerCell*meshHierarchy_(1).numElements);
	tmp = sort(edgeIndices,1)';
	[eNodMatTmp, ~, ~] = unique(tmp, 'rows');
	numEdgesTmp = size(eNodMatTmp,1);
	disp(['Number of Edges: ', sprintf('%d', numEdgesTmp)]);
end