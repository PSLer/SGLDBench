function SAGS_CallArora2019MatlabSuite_ExtractingGraph(resolution)
	global smoothedStressField_;
	global vertexEdgeGraph_;
	global frameStruct4Voxelization_;
	
	%% Extract the truss layout
	tStart = tic;
	truss = tex2CurvesTet(smoothedStressField_, round(resolution), false);
	disp(['Extracting Truss Layout Costs: ' sprintf('%10.3g',toc(tStart)) 's']);	
	% Generate sparser structure
	n = 1;  % generate sparse structure by selecting every n-th curve for each parameter (2==default)
	sparseTruss = sparsifyElements(truss, n);	
	% Generate manifold mesh using libigl's booleans
	% first, collapse tiny elements since the wire_mesh code does not tolerate
	% near-degenerate elements
	d = collapseSmallEdges(sparseTruss, 1e-3);
	node = d.Node;
	elem = d.Elem;
	% also, remove any disconnected nodes, possibly formed during the last operation
	[node, SVI] = remove_unreferenced(node, elem);
	elem = SVI(elem);

	vertexEdgeGraph_.state = 1;
	vertexEdgeGraph_.numNodes = size(node,1);
	vertexEdgeGraph_.nodeCoords = node;
	vertexEdgeGraph_.numEdges = size(elem,1);
	vertexEdgeGraph_.eNodMat = elem;	
	frameStruct4Voxelization_ = vertexEdgeGraph_;
	frameStruct4Voxelization_.edgeLengths = vecnorm(frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,1),:) ...
		- frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,2),:),2,2);		
end