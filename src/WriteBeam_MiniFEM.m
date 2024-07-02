function WriteBeam_MiniFEM(fileName, scalingEdgeThickness, nodeCoords_, edgeNodMat_)
	boundingBox_ = [min(nodeCoords_,[],1); max(nodeCoords_,[],1)];
	numNodes_ = size(nodeCoords_,1);
	numEles_ = size(edgeNodMat_,1);
	
	refDiameter = max(boundingBox_(2,:) - boundingBox_(1,:))/scalingEdgeThickness;
	diameterList_ = repmat(refDiameter, numEles_, 1);
	materialIndicatorField_ = ones(numEles_, 1);
	
	nodState_ = zeros(numNodes_,1);
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s ', 'Version');
	fprintf(fid, '%.1f\n', 2.0);
	
	fprintf(fid, '%s %s ', 'Frame Beam3D');
	fprintf(fid, '%d\n', 1);
	
	fprintf(fid, '%s ', 'Vertices:');
	fprintf(fid, '%d\n', numNodes_);		
	fprintf(fid, '%.6e %.6e %.6e\n', nodeCoords_');		
	fprintf(fid, '%s ', 'Elements:');
	fprintf(fid, '%d \n', numEles_);
	fprintf(fid, '%d %d %.6e %d\n', [double(edgeNodMat_) diameterList_ materialIndicatorField_]');	
	
	fprintf(fid, '%s %s ', 'Node State: ');
	fprintf(fid, '%d\n', numel(nodState_));
	if ~isempty(nodState_)
		fprintf(fid, '%d\n', nodState_);
	end

	fprintf(fid, '%s %s ', 'Node Forces:'); 
	fprintf(fid, '%d\n', 0);

	fprintf(fid, '%s %s ', 'Fixed Nodes:'); 
	fprintf(fid, '%d\n', 0);	
	fclose(fid);	
end