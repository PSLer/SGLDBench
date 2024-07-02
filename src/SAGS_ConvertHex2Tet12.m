function [hex_to_tet_nodeCoords, hex_to_tet_eNodMat, hex_to_tet_cartesianStress] = SAGS_ConvertHex2Tet12(solidMesh4Sim, cartesianStresses)
	numTetNodes = solidMesh4Sim.numNodes + solidMesh4Sim.numElements;
	eleCentX = solidMesh4Sim.nodeCoords(:,1); eleCentX = eleCentX(solidMesh4Sim.eNodMat);
	eleCentY = solidMesh4Sim.nodeCoords(:,2); eleCentY = eleCentY(solidMesh4Sim.eNodMat);
	eleCentZ = solidMesh4Sim.nodeCoords(:,3); eleCentZ = eleCentZ(solidMesh4Sim.eNodMat);
	eleCentroidList_ = [sum(eleCentX,2) sum(eleCentY,2) sum(eleCentZ,2)]/8;	
	
	hex_to_tet_nodeCoords = [solidMesh4Sim.nodeCoords; eleCentroidList_];
	numTetEles = 12*solidMesh4Sim.numElements;
	hex_to_tet_eNodMat = zeros(solidMesh4Sim.numElements,48);
	stressTensorAtElementCentroids = zeros(solidMesh4Sim.numElements, size(cartesianStresses,2));
	shapeFuncsAtCentroid = FEA_ShapeFunction(0.0, 0.0, 0.0);
	for ii=1:solidMesh4Sim.numElements
		iCartesianStress = cartesianStresses(solidMesh4Sim.eNodMat(ii,:), :);
		stressTensorAtElementCentroids(ii,:) = shapeFuncsAtCentroid * iCartesianStress;
	end
	hex_to_tet_cartesianStress = [cartesianStresses; stressTensorAtElementCentroids];
	
	newNodeIndices = (1:solidMesh4Sim.numElements)' + solidMesh4Sim.numNodes;
	hex_to_tet_eNodMat(:,1:4) = [solidMesh4Sim.eNodMat(:,[1 2 3]) newNodeIndices];
	hex_to_tet_eNodMat(:,5:8) = [solidMesh4Sim.eNodMat(:,[3 4 1]) newNodeIndices];
	hex_to_tet_eNodMat(:,9:12) = [solidMesh4Sim.eNodMat(:,[5 7 6]) newNodeIndices];
	hex_to_tet_eNodMat(:,13:16) = [solidMesh4Sim.eNodMat(:,[5 8 7]) newNodeIndices];
	
	hex_to_tet_eNodMat(:,17:20) = [solidMesh4Sim.eNodMat(:,[2 7 3]) newNodeIndices];
	hex_to_tet_eNodMat(:,21:24) = [solidMesh4Sim.eNodMat(:,[2 6 7]) newNodeIndices];
	hex_to_tet_eNodMat(:,25:28) = [solidMesh4Sim.eNodMat(:,[1 4 8]) newNodeIndices];
	hex_to_tet_eNodMat(:,29:32) = [solidMesh4Sim.eNodMat(:,[1 8 5]) newNodeIndices];
	
	hex_to_tet_eNodMat(:,33:36) = [solidMesh4Sim.eNodMat(:,[1 5 2]) newNodeIndices];			
	hex_to_tet_eNodMat(:,37:40) = [solidMesh4Sim.eNodMat(:,[5 6 2]) newNodeIndices];
	hex_to_tet_eNodMat(:,41:44) = [solidMesh4Sim.eNodMat(:,[8 4 3]) newNodeIndices];
	hex_to_tet_eNodMat(:,45:48) = [solidMesh4Sim.eNodMat(:,[7 8 3]) newNodeIndices];
	hex_to_tet_eNodMat = hex_to_tet_eNodMat'; hex_to_tet_eNodMat = hex_to_tet_eNodMat(:);
	hex_to_tet_eNodMat = reshape(hex_to_tet_eNodMat, 4, numTetEles)';
end