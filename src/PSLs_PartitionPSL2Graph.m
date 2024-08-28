function [iGraphNodes, iGraphEdges] = PSLs_PartitionPSL2Graph(PSLcoords, graphSpan, numExisingNode)
	numIntPoints = size(PSLcoords,1);
	iGraphNodeIndices = 1:graphSpan:numIntPoints;
	if numIntPoints-iGraphNodeIndices(end)>1
		iGraphNodeIndices(end+1) = numIntPoints;
	end
	numiGraphNodes = numel(iGraphNodeIndices);    
	iGraphNodes = PSLcoords(iGraphNodeIndices,:);
	iGraphNodeIndicesLocal = 1:numel(iGraphNodeIndices);
	if numiGraphNodes > 2
		iGraphEdges = [iGraphNodeIndicesLocal(1:end-1); iGraphNodeIndicesLocal(2:end)]';
	else
		iGraphEdges = [1 2];
	end
	iGraphEdges = iGraphEdges + numExisingNode;
end