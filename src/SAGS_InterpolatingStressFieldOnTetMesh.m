function SAGS_InterpolatingStressFieldOnTetMesh()
	global outPath_;
	global meshHierarchy_;
	global nodeCoords_;
	global cartesianStressField_;	
	global stressFieldOnTetMesh_; 
	global gateWayTetMesh_;
	global dataPrep4SAGS_;
	
	nodeCoords_ = zeros(meshHierarchy_(1).numNodes,3);
	nodeCoords_(:,1) = double(niftiread(strcat(outPath_, 'cache_coordX.nii')));
	nodeCoords_(:,2) = double(niftiread(strcat(outPath_, 'cache_coordY.nii')));
	nodeCoords_(:,3) = double(niftiread(strcat(outPath_, 'cache_coordZ.nii')));
	
	stressFieldOnTetMesh_ = zeros(gateWayTetMesh_.numNodes,6);
	for ii=1:gateWayTetMesh_.numNodes
		iNode = gateWayTetMesh_.nodeCoords(ii,:);
		[eleIndex, paraCoordinates, opt] = Common_LocatePointOnCartesianMesh(iNode);
		if opt
			NIdx = meshHierarchy_(1).eNodMat(eleIndex,:);
			eleCartesianStress = cartesianStressField_(NIdx,:);
			stressFieldOnTetMesh_(ii,:) = ...
				FEA_ShapeFunction(paraCoordinates(1), paraCoordinates(2), paraCoordinates(3)) * eleCartesianStress;					
		else
			[~, closestBoundaryNode] = min(vecnorm(iNode-meshHierarchy_(1).boundaryNodeCoords,2,2));
			closestNode =  meshHierarchy_(1).nodesOnBoundary(closestBoundaryNode);
			stressFieldOnTetMesh_(ii,:) = cartesianStressField_(closestNode,:);
		end
	end

	dataPrep4SAGS_ = struct('nodeCoords', [], 'eNodMat', [], 'cartesianStress', [], 'ps', [], 'vM', [], 'frameField', []);
	dataPrep4SAGS_.nodeCoords = gateWayTetMesh_.nodeCoords;
	dataPrep4SAGS_.eNodMat = gateWayTetMesh_.eNodMat;
	dataPrep4SAGS_.cartesianStress = stressFieldOnTetMesh_;
	dataPrep4SAGS_.vM = sqrt(0.5*((stressFieldOnTetMesh_(:,1)-stressFieldOnTetMesh_(:,2)).^2 + ...
		(stressFieldOnTetMesh_(:,2)-stressFieldOnTetMesh_(:,3)).^2 + (stressFieldOnTetMesh_(:,3)...
			-stressFieldOnTetMesh_(:,1)).^2 ) + 3*( stressFieldOnTetMesh_(:,6).^2 + stressFieldOnTetMesh_(:,4).^2 + ...
				stressFieldOnTetMesh_(:,5).^2 ));	
	dataPrep4SAGS_.ps = zeros(size(dataPrep4SAGS_.nodeCoords,1),12);
	for ii=1:size(dataPrep4SAGS_.nodeCoords,1)
		dataPrep4SAGS_.ps(ii,:)= FEA_ComputePrincipalStress(dataPrep4SAGS_.cartesianStress(ii,:));
	end	
	
	%%Write Gateway Stress File for Graded Voronoi Diagram Generation
	IO_ExportStressField2TSV();
end

function [nextElementIndex, paraCoordinates, opt] = Common_LocatePointOnCartesianMesh(physicalCoordinates)
	global meshHierarchy_;
	global nodeCoords_;
	global boundingBox_;
	
	nextElementIndex = 0; paraCoordinates = []; opt = 0;
	resX = meshHierarchy_(1).resX;
	resY = meshHierarchy_(1).resY;
	resZ = meshHierarchy_(1).resZ;	
	physicalCoordinates = physicalCoordinates - boundingBox_(1,:);
	if 0==physicalCoordinates(1)
		eleX = 1;				
	else
		eleX = ceil(physicalCoordinates(1)/meshHierarchy_(1).eleSize(1));
		if eleX<1 || eleX>resX, return; end
	end
	if 0==physicalCoordinates(2)
		eleY = 1;
	else
		eleY = ceil(physicalCoordinates(2)/meshHierarchy_(1).eleSize(2));
		if eleY<1 || eleY>resY, return; end
	end
	if 0==physicalCoordinates(3)
		eleZ = 1;
	else
		eleZ = ceil(physicalCoordinates(3)/meshHierarchy_(1).eleSize(3));
		if eleZ<1 || eleZ>resZ, return; end
	end			
	
	tarEle = resX*resY*(eleZ-1) + resY*(eleX-1)+(resY-eleY+1);
	nextElementIndex = meshHierarchy_(1).eleMapForward(tarEle);
	if nextElementIndex	
		opt = 1;
		% relatedNodes = Common_RecoverHalfeNodMat(meshHierarchy_(1).eNodMatHalf(nextElementIndex,:));
		relatedNodes = meshHierarchy_(1).eNodMat(nextElementIndex,:);
		relatedNodeCoords = nodeCoords_(relatedNodes',:)-boundingBox_(1,:);
		paraCoordinates = 2*(physicalCoordinates - relatedNodeCoords(1,:)) / meshHierarchy_(1).eleSize(1) - 1;
	end	
end