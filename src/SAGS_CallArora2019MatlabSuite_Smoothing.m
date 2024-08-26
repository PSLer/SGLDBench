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
end