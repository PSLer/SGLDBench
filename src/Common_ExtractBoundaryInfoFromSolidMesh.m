function [reOrderedBoundaryNodeCoords, reOrderedBoundaryFaceNodMat, nodState, eleState] = ...
			Common_ExtractBoundaryInfoFromSolidMesh(nodeCoords, eNodMat)
	
	numNodes = size(nodeCoords,1);
	[numEles, numNodesPerEle] = size(eNodMat);
	switch numNodesPerEle
		case 4
			numNodesPerFace = 3;
			numFacesPerEle = 4;
			numEdgesPerEle = 6;
			patchIndices = eNodMat(:, [3 2 1  1 2 4  2 3 4  3 1 4])';
		case 8
			numNodesPerFace = 4;
			numFacesPerEle = 6;
			numEdgesPerEle = 12;
			patchIndices = eNodMat(:, [4 3 2 1  5 6 7 8  1 2 6 5  8 7 3 4  5 8 4 1  2 3 7 6])';
	end
	patchIndices = reshape(patchIndices(:), numNodesPerFace, numFacesPerEle*numEles)';	
	tmp = sort(patchIndices,2);
	[uniqueFaces, ia, ic] = unique(tmp, 'stable', 'rows');
	leftFaceIDs = (1:numFacesPerEle*numEles)'; leftFaceIDs = setdiff(leftFaceIDs, ia);
	leftFaces = tmp(leftFaceIDs,:);
	[surfFaces, surfFacesIDsInUniqueFaces] = setdiff(uniqueFaces, leftFaces, 'rows');
	boundaryFaceNodMat = patchIndices(ia(surfFacesIDsInUniqueFaces),:);
	boundaryNodes = unique(boundaryFaceNodMat);
	nodState = zeros(numNodes,1); nodState(boundaryNodes) = 1;	
	eleState = numEdgesPerEle*ones(numEles,1);	
	
	allNodes = zeros(numNodes, 1);
	numBoundaryNodes = numel(boundaryNodes);
	allNodes(boundaryNodes,1) = (1:numBoundaryNodes)';
	reOrderedBoundaryFaceNodMat = allNodes(boundaryFaceNodMat);
	reOrderedBoundaryNodeCoords = nodeCoords(boundaryNodes,:);
end