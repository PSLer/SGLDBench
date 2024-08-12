function [passiveElementsOnBoundary, passiveElementsNearLoads, passiveElementsNearFixation] = TopOpti_SetPassiveElements(numLayerboundary, numLayerLoads, numLayerFixation)
	global meshHierarchy_;
	global loadingCond_;
	global fixingCond_;
	global passiveElements_; 
	
	numLayerboundary = round(numLayerboundary);
	numLayerLoads = round(numLayerLoads);
	numLayerFixation = round(numLayerFixation);
	existingPassiveElements = passiveElements_;
    
	passiveElementsOnBoundary = [];
	if numLayerboundary>0
		index = 1;
		while index<=numLayerboundary
			if 1==index
				passiveElementsOnBoundary = double(meshHierarchy_(1).elementsOnBoundary);
			else
				passiveElementsOnBoundary = Common_IncludeAdjacentElements(passiveElementsOnBoundary);
			end
			index = index + 1;
		end		
	end
	passiveElementsNearLoads = [];
	passiveElementsNearFixation = [];
	if numLayerLoads>0 || numLayerFixation>0
		%%Relate Elements to Nodes
		allElements = zeros(meshHierarchy_(1).numElements,1);
		allElements(meshHierarchy_(1).elementsOnBoundary) = 1;
		nodeStruct_ = struct('arr', []);
		nodeStruct_ = repmat(nodeStruct_, meshHierarchy_(1).numNodes, 1);
		boundaryNodes_Temp = [];
		numNodsPerEle = 8;
		for ii=1:meshHierarchy_(1).numElements
			if 0 || allElements(ii) %%switch 0 to 1 for non-boundary fixation situations, efficiency loss
				% iNodes = meshHierarchy_(1).eNodMatHalf(ii,:);
				% iNodes = Common_RecoverHalfeNodMat(iNodes);
				iNodes = meshHierarchy_(1).eNodMat(ii,:);
				for jj=1:numNodsPerEle
					nodeStruct_(iNodes(jj)).arr(1,end+1) = ii;
				end
			end
		end
		
		%% Extract Elements near Loads
		passiveElementsNearLoads = [];
		if numLayerLoads>0
			loadedNodes = meshHierarchy_(1).nodesOnBoundary(loadingCond_(:,1));
			allLoadedNodes = nodeStruct_(loadedNodes);
			passiveElementsNearLoads = unique([allLoadedNodes.arr]);
			index = 2;
			while index<=numLayerLoads
				passiveElementsNearLoads = Common_IncludeAdjacentElements(passiveElementsNearLoads);
				index = index + 1;
			end				
		end
		passiveElementsNearFixation = [];
		if numLayerFixation>0
			fixedNodes = meshHierarchy_(1).nodesOnBoundary(fixingCond_(:,1));
			allFixedNodes = nodeStruct_(fixedNodes);
			passiveElementsNearFixation = unique([allFixedNodes.arr]);
			index = 2;
			while index<=numLayerFixation
				passiveElementsNearFixation = Common_IncludeAdjacentElements(passiveElementsNearFixation);
				index = index + 1;
			end
		end
		
		
		passiveElements_ = unique([existingPassiveElements; passiveElementsOnBoundary(:); ...
			passiveElementsNearLoads(:); passiveElementsNearFixation(:)]);
	else
		passiveElements_ = unique([existingPassiveElements; passiveElementsOnBoundary(:)]);
	end
	passiveElements_ = unique([existingPassiveElements; passiveElements_]);
end