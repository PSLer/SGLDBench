function SAGS_CallGao2017ExecutableProgram(aspectRatio)
	global dataPrep4SAGS_;
	global vertexEdgeGraph_;
	global frameStruct4Voxelization_;
	
	%% Initialize Frame Field
	aspectRatio = max([1.0, aspectRatio]);
	dataPrep4SAGS_.frameField = dataPrep4SAGS_.ps(:, [2 3 4 6 7 8 10 11 12 1 5 9]);
	if 1==aspectRatio
		dataPrep4SAGS_.frameField(:,[10 11 12]) = ones(size(dataPrep4SAGS_.frameField,1),3);
	else
		minVal = 1/aspectRatio;
		dataPrep4SAGS_.frameField(:,[10 11 12]) = abs(dataPrep4SAGS_.frameField(:,[10 11 12]));
		for ii=1:size(dataPrep4SAGS_.nodeCoords,1)
			iFrame = abs(dataPrep4SAGS_.frameField(ii,[10 11 12]));
			iFrame = iFrame / max(iFrame);
			iFrame = max(iFrame, minVal);
			dataPrep4SAGS_.frameField(ii,[10 11 12]) = iFrame;
		end	
	end
	dataPrep4SAGS_.frameField = dataPrep4SAGS_.frameField';
	dataPrep4SAGS_.frameField = dataPrep4SAGS_.frameField(:);
	dataPrep4SAGS_.frameField = reshape(dataPrep4SAGS_.frameField, 3, 4*size(dataPrep4SAGS_.nodeCoords,1))';
	
	%% Write Data into 
	%% Mesh
	fid = fopen('./out/FrameData4Gao2017.mesh', 'w');	
	fprintf(fid, '%d %d\n', [size(dataPrep4SAGS_.nodeCoords,1) size(dataPrep4SAGS_.eNodMat,1)]);
	fprintf(fid, '%.6f %.6f %.6f\n', dataPrep4SAGS_.nodeCoords');
	fprintf(fid, '%d %d %d %d %d\n', [4*ones(size(dataPrep4SAGS_.eNodMat,1),1) dataPrep4SAGS_.eNodMat-1]');
	fclose(fid);
	%%Frame Field
	fid = fopen('./out/FrameData4Gao2017.txt', 'w');	
	fprintf(fid, '%.6f %.6f %.6f\n', dataPrep4SAGS_.frameField');
	fclose(fid);	
	
	%%Run Gao2017
	system('"./src/external/Gao2017/tensor-field-meshing.exe" -b -i ./out/FrameData4Gao2017');
	
	%%Read Field-aligned Graph in
	IO_ImportVertexEdgeGraph('./out/FrameData4Gao2017_graph_opt.obj');
	
	frameStruct4Voxelization_ = vertexEdgeGraph_;
	frameStruct4Voxelization_.edgeLengths = vecnorm(frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,1),:) ...
		- frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,2),:),2,2);		
end