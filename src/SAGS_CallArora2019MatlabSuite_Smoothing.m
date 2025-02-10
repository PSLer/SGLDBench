function SAGS_CallArora2019MatlabSuite_Smoothing()
	global dataPrep4SAGS_;
	global smoothedStressField_;
	
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
	
	smoothedStressField_ = dataOut;
	
	IO_ExportStressField2TSV_eleWise();
end

function IO_ExportStressField2TSV_eleWise()
	global outPath_;
	global smoothedStressField_; 
	global gateWayTetMesh_;	
	
	fileName = strcat(outPath_, 'StressField_Tet_v2_eleWise.stress');
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s ', 'Version');
	fprintf(fid, '%.1f\n', 2.0);
	
	fprintf(fid, '%s %s ', 'Solid Tet');
	fprintf(fid, '%d\n', 1);
	
	fprintf(fid, '%s ', 'Vertices:');
	fprintf(fid, '%d\n', gateWayTetMesh_.numNodes);		
	fprintf(fid, '%.6e %.6e %.6e\n', gateWayTetMesh_.nodeCoords');

	fprintf(fid, '%s ', 'Elements:');
	fprintf(fid, '%d \n', gateWayTetMesh_.numElements);
	fprintf(fid, '%d %d %d %d\n', gateWayTetMesh_.eNodMat');

	fprintf(fid, '%s %s ', 'Node Forces:'); 
	fprintf(fid, '%d\n', 0);
	fprintf(fid, '%s %s ', 'Fixed Nodes:'); fprintf(fid, '%d\n', 0);

	fprintf(fid, '%s %s', 'Cartesian Stress:'); 
	fprintf(fid, '%d\n', gateWayTetMesh_.numElements);
	dirField = ones(gateWayTetMesh_.numElements, 12);
	dirField(:,2:4) = smoothedStressField_.t;
	dirField(:,6:8) = smoothedStressField_.w;
	dirField(:,10:12) = smoothedStressField_.v;
	fprintf(fid, '%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n', dirField');
	fclose(fid);
end