function SAGS_CallArora2019MatlabSuite(resolution)
	global dataPrep4SAGS_;
	global vertexEdgeGraph_;
	global frameStruct4Voxelization_;
	
	%% Initialize Stress Field
	inputData = struct('V', [], 'T', [], 'stress', []);
	inputData.V = dataPrep4SAGS_.nodeCoords;
	inputData.T = dataPrep4SAGS_.eNodMat;
	stressField = dataPrep4SAGS_.cartesianStress(:,[1 2 3 6 4 5]);
	shapeFuncsAtCentroid = [0.25 0.25 0.25 0.25];
	inputData.stress = zeros(size(inputData.T,1),6);
	for ii=1:size(inputData.T)
		iStress = stressField(inputData.T(ii,:),:);
		inputData.stress(ii,:) = shapeFuncsAtCentroid * iStress;
	end		
	
	%%Smoothing
	tStart = tic;
	dataFrames = fitFramesToData3D(inputData.V, inputData.T, inputData.stress);
	disp(['Smoothing Stress Field Costs: ' sprintf('%10.3g',toc(tStart)) 's']);	
	%% Solve for texture parametrization
	tStart = tic;
	b = 0.1;  % select `beta=0.1` for parametrization solve
	dataTex = fitTexCoords3D(dataFrames, b);
	dataOut = dataTex;
	dataOut.u = matrixnormalize(dataOut.u);
	disp(['Solving Texture Parametrization Costs: ' sprintf('%10.3g',toc(tStart)) 's']);	
	%% Extract the truss layout
	tStart = tic;
	truss = tex2CurvesTet(dataOut, resolution, false);
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