function TopOpti_BuildDensityFilter_matrixFree()
	global outPath_;
	global meshHierarchy_;
	global uniqueCellsInDensityFilteringMapVec_;
	global adjInfoUniqueCellDensityFiltering_;
	global sumWeightsDensityFilter_;
	global identicalWeightsDensityFiltering_;
	global identicalCellAdjInfo_;
	global rMin_;

	%%1. identify unique cells in filtering (boundary cells)
	numElements = meshHierarchy_(1).numElements;
	numLayerboundary = ceil(rMin_);
	index = 1;
	while index<=numLayerboundary
		if 1==index
			uniqueCellsInDensityFiltering = double(meshHierarchy_(1).elementsOnBoundary);
		else
			uniqueCellsInDensityFiltering = Common_IncludeAdjacentElements(uniqueCellsInDensityFiltering);
		end
		index = index + 1;
	end
	uniqueCellsInDensityFiltering = unique(uniqueCellsInDensityFiltering);
	numUniqueCells = numel(uniqueCellsInDensityFiltering);
	uniqueCellsInDensityFilteringMapVec_ = zeros(meshHierarchy_(1).numElements, 1, 'int32');
	uniqueCellsInDensityFilteringMapVec_(uniqueCellsInDensityFiltering) = (1:numUniqueCells)';
	% uniqueWeightsSum = zeros(numUniqueCells,1,'single');
	uniqueWeightsSum = zeros(numUniqueCells,1);
	adjInfoUniqueCellDensityFiltering_ = struct('cells', [], 'weights', []);
	adjInfoUniqueCellDensityFiltering_ = repmat(adjInfoUniqueCellDensityFiltering_, numUniqueCells, 1);
	
	%%2. Build Adj-info for unique cells
	%%	1	4	7		10	 13	  16		19	 22	  25
	%%	2	5	8		11	 14*  17		20	 23	  26
	%%	3	6	9		12	 15   18		21	 24	  27
	%%	 bottom				middle				top
	resX = meshHierarchy_(1).resX;
	resY = meshHierarchy_(1).resY;
	resZ = meshHierarchy_(1).resZ;
	eleMapForward = meshHierarchy_(1).eleMapForward;
	eleMapBack = meshHierarchy_(1).eleMapBack;
	% eleSize = single(meshHierarchy_(1).eleSize(1));
	eleSize = meshHierarchy_(1).eleSize(1);
	
	% eleCentroidList = niftiread(strcat(outPath_, 'cache_eleCentroidList.nii'));
	% eleCentroidList = double(eleCentroidList);
	rMin = rMin_;
	zeroBolck = zeros((2*ceil(rMin)-1)^3,1);
	for kk = 1:resZ
		for ii = 1:resX
			for jj = 1:resY			
				e1MapBack = (kk-1)*resX*resY+(ii-1)*resY+jj;
				e1 = eleMapForward(e1MapBack);
                if ~e1, continue; end
				e1Unique = uniqueCellsInDensityFilteringMapVec_(e1);
				if e1Unique					
					adjCells = zeroBolck;
					e1Weights = zeroBolck;
					iIndex = 0;
					for kk2 = max(kk-(ceil(rMin)-1),1):min(kk+(ceil(rMin)-1),resZ)
						for ii2 = max(ii-(ceil(rMin)-1),1):min(ii+(ceil(rMin)-1),resX)
							for jj2 = max(jj-(ceil(rMin)-1),1):min(jj+(ceil(rMin)-1),resY)
								e2MapBack = (kk2-1)*resX*resY+(ii2-1)*resY+jj2;
								e2 = eleMapForward(e2MapBack);
								if e2
									iIndex = iIndex+1;
									adjCells(iIndex,1) = e2;
									e1Weights(iIndex,1) = rMin - norm([ii jj kk]-[ii2 jj2 kk2]);
								end
							end
						end
					end
					adjCells = adjCells(1:iIndex,1);
					e1Weights = e1Weights(1:iIndex,1);
					% e1Weights = rMin - norm([ii jj kk]-[ii2 jj2 kk2]);
					% e1Weights = rMin*eleSize - vecnorm(eleCentroidList(e1,:)-eleCentroidList(adjCells,:),2,2);
					e1Weights(e1Weights<0) = 0;
					adjInfoUniqueCellDensityFiltering_(e1Unique).cells = adjCells;
					adjInfoUniqueCellDensityFiltering_(e1Unique).weights = e1Weights(:);
					uniqueWeightsSum(e1Unique) = sum(e1Weights);
				end
			end
		end
	end

	for ii=1:numUniqueCells
		iCellsMapForward = adjInfoUniqueCellDensityFiltering_(ii).cells;
		iCellsMapBack = eleMapBack(iCellsMapForward);
		adjInfoUniqueCellDensityFiltering_(ii).cells = iCellsMapBack;
	end

	%%3. Build Adj-info for cells with identical adj-topology	
	cond = 0;
	for kk = 1:resZ
		for ii = 1:resX
			for jj = 1:resY	
				expIdenticalCellMapBack = (kk-1)*resX*resY+(ii-1)*resY+jj;
				expIdenticalCell = eleMapForward(expIdenticalCellMapBack);				
				if ~expIdenticalCell, continue; end
				e1Identical = uniqueCellsInDensityFilteringMapVec_(expIdenticalCell);
				if 0==e1Identical				
					adjCellsIdentical = zeroBolck;
					expIdenticalCellWeights = zeroBolck;
					iIndex = 0;
					for kk2 = max(kk-(ceil(rMin)-1),1):min(kk+(ceil(rMin)-1),resZ)
						for ii2 = max(ii-(ceil(rMin)-1),1):min(ii+(ceil(rMin)-1),resX)
							for jj2 = max(jj-(ceil(rMin)-1),1):min(jj+(ceil(rMin)-1),resY)
								e2MapBack = (kk2-1)*resX*resY+(ii2-1)*resY+jj2;
								e2 = eleMapForward(e2MapBack);
								if e2
									iIndex = iIndex+1;
									adjCellsIdentical(iIndex,1) = e2;
									expIdenticalCellWeights(iIndex,1) = rMin - norm([ii jj kk]-[ii2 jj2 kk2]);
								end
							end
						end
					end					
					adjCellsIdentical = adjCellsIdentical(1:iIndex,1);
					expIdenticalCellWeights = expIdenticalCellWeights(1:iIndex,1);
					% expIdenticalCellWeights = rMin*eleSize - vecnorm(eleCentroidList(expIdenticalCell,:)-eleCentroidList(adjCellsIdentical,:),2,2);		
					expIdenticalCellWeights(expIdenticalCellWeights<0) = 0;
					cond = 1; break;
				end				
			end
			if cond, break; end
		end
		if cond, break; end
	end
	identicalWeightsDensityFiltering_ = expIdenticalCellWeights(:);
	identicalCellAdjInfo_ = eleMapBack(adjCellsIdentical) - expIdenticalCellMapBack;
	% identicalCellAdjInfo_ = identicalCellAdjInfo_(:)';
	
	%identicalWeightsSum = single(sum(expIdenticalCellWeights));
	identicalWeightsSum = sum(expIdenticalCellWeights);
	sumWeightsDensityFilter_ = repmat(identicalWeightsSum, numElements, 1);
	sumWeightsDensityFilter_(uniqueCellsInDensityFiltering) = uniqueWeightsSum;
end