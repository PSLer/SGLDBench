function IO_ImportVertexEdgeGraph(fileName)
	[~,~,dataType] = fileparts(fileName);
	switch dataType
		case '.xxx'
			%%...
		case '.obj'
			IO_ImportVertexEdgeGraph_Format_obj(fileName);
	end
end

function IO_ImportVertexEdgeGraph_Format_obj(fileName)
	global vertexEdgeGraph_;
	%%Read Data
	nodeCoords = []; eNodMat = [];
	fid = fopen(fileName, 'r');
	while 1
		tline = fgetl(fid);
		if ~ischar(tline),   break,   end  % exit at end of file 
		ln = sscanf(tline,'%s',1); % line type 
		switch ln
			case 'v' % graph vertexs
				nodeCoords(end+1,1:3) = sscanf(tline(2:end), '%f')';
			case 'l'
				eNodMat(end+1,1:2) = sscanf(tline(2:end), '%d')';
		end
	end
	fclose(fid);
	vertexEdgeGraph_.nodeCoords = nodeCoords;
	vertexEdgeGraph_.eNodMat = eNodMat;
	vertexEdgeGraph_.numNodes = size(vertexEdgeGraph_.nodeCoords,1); 
	vertexEdgeGraph_.numEdges = size(vertexEdgeGraph_.eNodMat,1);
	vertexEdgeGraph_.edgeLengths = vecnorm(nodeCoords(eNodMat(:,1),:) - nodeCoords(eNodMat(:,2),:),2,2);
	vertexEdgeGraph_.state = 1;	
end