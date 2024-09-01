function aveDensVec = TopOpti_DensityFiltering_matrixFree(densVec, opt)
	global meshHierarchy_;
	global voxelizedVolume_;
	global uniqueCellsInDensityFilteringMapVec_;
	global adjInfoUniqueCellDensityFiltering_;	
	global identicalWeightsDensityFiltering_;
	global identicalCellAdjInfo_;
	global sumWeightsDensityFilter_;
	
	aveDensVec = densVec;
	uniqueCellsInDensityFilteringMapVec = uniqueCellsInDensityFilteringMapVec_;
	adjInfoUniqueCellDensityFiltering = adjInfoUniqueCellDensityFiltering_;
	identicalCellAdjInfo = identicalCellAdjInfo_;
	identicalWeightsDensityFiltering = identicalWeightsDensityFiltering_;
	sumWeightsDensityFilter = sumWeightsDensityFilter_;
	numElements = meshHierarchy_(1).numElements;
	eleMapBack = meshHierarchy_(1).eleMapBack;
	% rhoMap = zeros(size(meshHierarchy_(1).eleMapForward), 'single');
	rhoMap = zeros(size(meshHierarchy_(1).eleMapForward));
	switch opt
		case 1
			% aveDensVec = H_*(densVec./Hs_);
			densVec = densVec ./ sumWeightsDensityFilter;
			rhoMap(voxelizedVolume_) = densVec;		
			if isempty(gcp('nocreate')), parpool('Threads'); end		
			parfor ii=1:numElements
				iUniqueEle = uniqueCellsInDensityFilteringMapVec(ii);
				if iUniqueEle
					adjCellInfo = adjInfoUniqueCellDensityFiltering(iUniqueEle);
					adjCellDens = rhoMap(adjCellInfo.cells);
					aveDensVec(ii) = adjCellDens(:)' * adjCellInfo.weights;
				else
					iEleMapBack = eleMapBack(ii);
					adjCellsMapBack = identicalCellAdjInfo + iEleMapBack;
					adjCellDens = rhoMap(adjCellsMapBack);
					aveDensVec(ii) = adjCellDens(:)' * identicalWeightsDensityFiltering;
				end
			end
		case 0
			% aveDensVec = H_*densVec./Hs_;
			rhoMap(voxelizedVolume_) = densVec;			
			if isempty(gcp('nocreate')), parpool('Threads'); end					
			parfor ii=1:numElements
				iUniqueEle = uniqueCellsInDensityFilteringMapVec(ii);
				if iUniqueEle
					adjCellInfo = adjInfoUniqueCellDensityFiltering(iUniqueEle);				
					adjCellDens = rhoMap(adjCellInfo.cells);
					aveDensVec(ii) = adjCellDens(:)' * adjCellInfo.weights;
				else
					iEleMapBack = eleMapBack(ii);
					adjCellsMapBack = identicalCellAdjInfo + iEleMapBack;
					adjCellDens = rhoMap(adjCellsMapBack);
					aveDensVec(ii) = adjCellDens(:)' * identicalWeightsDensityFiltering;
				end
			end
			aveDensVec = aveDensVec ./ sumWeightsDensityFilter;
		otherwise
			error('Wrong option for checker board filtering')
	end	
end