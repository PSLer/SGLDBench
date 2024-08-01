function PSLs_GeneratePSLsBy3DTSV(lineDensCtrl, psDirIndicator)
	global boundingBox_;
	global meshHierarchy_;

	global mergeTrigger_;
	global selectedPrincipalStressField_;
	global tracingStepWidth_;
	global integrationStepLimit_;
	global relaxedFactor_;
	global permittedMaxAdjacentTangentAngleDeviation_;
	global minLengthVisiblePSLs_;
	
	global eleSize_;	
	
	%%To facilitate PSLs tracing
	eleSize_ = min(meshHierarchy_(1).eleSize);
	mergeTrigger_ = min(boundingBox_(2,:)-boundingBox_(1,:))/lineDensCtrl;
	selectedPrincipalStressField_ = psDirIndicator;
	tracingStepWidth_ = 0.5 * eleSize_;
	integrationStepLimit_ = ceil(1.5*norm(boundingBox_(2,:)-boundingBox_(1,:))/tracingStepWidth_);
	permittedMaxAdjacentTangentAngleDeviation_ = 20; %% pi/permittedMaxAdjacentTangentAngleDeviation_
	relaxedFactor_ = 1.0;
	minLengthVisiblePSLs_ = 20;
	
	GenerateSpaceFillingPSLs();
end

function GenerateSpaceFillingPSLs()
	global outPath_;
	global meshHierarchy_;
	global mergeTrigger_;
	global seedPoints_;
    global seedPointsHistory_;
	global seedPointsValence_; 
	global selectedPrincipalStressField_;
	global majorPSLpool_; 
    global mediumPSLpool_; 
    global minorPSLpool_; 
	global majorCoordList_; 
    global mediumCoordList_; 
    global minorCoordList_;
    global relaxedFactor_; 
	
	%%Initialize Seed Points
	seedDensCtrl = ceil(mergeTrigger_/sqrt(3) * 1.0); seedDensCtrl = max(seedDensCtrl, 4);
	validElementVolume = reshape(meshHierarchy_(1).eleMapForward, meshHierarchy_(1).resY, ...
		meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	
	sampledElements = validElementVolume(seedDensCtrl:seedDensCtrl:meshHierarchy_(1).resY-seedDensCtrl, ...
		seedDensCtrl:seedDensCtrl:meshHierarchy_(1).resX-seedDensCtrl, seedDensCtrl:seedDensCtrl:meshHierarchy_(1).resZ-seedDensCtrl);	
	sampledElements = sampledElements(:);
	sampledElements(0==sampledElements) = [];
	% seedPointsHistory_ = meshHierarchy_(1).eleCentroidList(sampledElements,:);
	eleCentroidList = niftiread(strcat(outPath_, 'cache_eleCentroidList.nii'));
	seedPointsHistory_ = double(eleCentroidList(sampledElements,:));
	startCoord_ = (min(seedPointsHistory_, [], 1) + max(seedPointsHistory_, [], 1)) / 2;
    seedPoints_ = seedPointsHistory_;
	numSeedPoints = size(seedPoints_,1);	
    seedPointsValence_ = ones(numSeedPoints, 3);
	if selectedPrincipalStressField_(1), seedPointsValence_(:,1) = 0; end
	if selectedPrincipalStressField_(2), seedPointsValence_(:,2) = 0; end
	if selectedPrincipalStressField_(3), seedPointsValence_(:,3) = 0; end
	numPSF = numel(1==selectedPrincipalStressField_);
	majorPSLpool_ = Data_PrincipalStressLineStruct();
	mediumPSLpool_ = Data_PrincipalStressLineStruct();
	minorPSLpool_ = Data_PrincipalStressLineStruct();	
	
	%% Iteration
	its = 0;
	looper = sum(sum(seedPointsValence_));	
	while looper<3*numSeedPoints
		its = its + 1;
		valenceMetric = sum(seedPointsValence_,2);
		%% 1st Priority: semi-empty seeds > empty seeds, which helps get PSLs intersection
		%% 2nd Priority: seeds with same valence, the one closest to the start point goes first
		switch numPSF
			case 1
				unFinishedSppsValence2 = find(2==valenceMetric);
				[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence2,end-2:end),2,2));
				spp = unFinishedSppsValence2(tarPos);
			case 2
				unFinishedSppsValence2 = find(2==valenceMetric); 
				if ~isempty(unFinishedSppsValence2) %% 1st Priority
					[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence2,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence2(tarPos);
				else
					unFinishedSppsValence1 = find(1==valenceMetric);
					[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence1,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence1(tarPos);			
				end					
			case 3
				unFinishedSppsValence12 = find(3>valenceMetric); 
				unFinishedSppsValence12 = unFinishedSppsValence12(valenceMetric(unFinishedSppsValence12)>0); 
				if ~isempty(unFinishedSppsValence12) %% 1st Priority
					[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence12,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence12(tarPos);
				else
					unFinishedSppsValence0 = find(0==valenceMetric);
					[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence0,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence0(tarPos);		
				end						
		end
		valences = seedPointsValence_(spp,:);						
		seed = seedPoints_(spp,:);
		if 0==valences(1)
			seedPointsValence_(spp,1) = 1;
			majorPSL = CreatePrincipalStressLine(seed, 'MAJOR');			
			if 0==majorPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end			
			majorPSLpool_(end+1,1) = majorPSL;				
			majorCoordList_(end+1:end+majorPSL.length,:) = majorPSL.phyCoordList;
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
			if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve(seedPoints_(...
						sppsEmptyMajorValence,:), majorPSL.phyCoordList);					
				potentialSolidSppsMajor = find(potentialDisListMajor<relaxedFactor_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);								
					seedPointsValence_(spps2BeMerged,1) = 1;
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValence_(modifiedMediumValences,2) = 1;					
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');				
					seedPointsValence_(modifiedMinorValences,3) = 1;	
				end
			end				
		end
		
		if 0==valences(2)
			seedPointsValence_(spp,2) = 1;
			mediumPSL = CreatePrincipalStressLine(seed, 'MEDIUM');
			if 0==mediumPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end
			mediumPSLpool_(end+1,1) = mediumPSL;
			mediumCoordList_(end+1:end+mediumPSL.length,:) = mediumPSL.phyCoordList;
			sppsEmptyMediumValence = find(0==seedPointsValence_(:,2));
			if ~isempty(sppsEmptyMediumValence)
				[potentialDisListMedium, potentialPosListMedium] = GetDisListOfPointList2Curve(seedPoints_(...
						sppsEmptyMediumValence,:), mediumPSL.phyCoordList);					
				potentialSolidSppsMedium = find(potentialDisListMedium<relaxedFactor_);
				if ~isempty(potentialSolidSppsMedium)
					spps2BeMerged = sppsEmptyMediumValence(potentialSolidSppsMedium);
					seedPoints_(spps2BeMerged,:) = potentialPosListMedium(potentialSolidSppsMedium,:);								
					seedPointsValence_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValence_(modifiedMajorValences,1) = 1;
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValence_(modifiedMinorValences,3) = 1;
				end
			end				
		end		
		
		if 0==valences(3)
			seedPointsValence_(spp,3) = 1;			
			minorPSL = CreatePrincipalStressLine(seed, 'MINOR');
			if 0==minorPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end		
			minorPSLpool_(end+1,1) = minorPSL;
			minorCoordList_(end+1:end+minorPSL.length,:) = minorPSL.phyCoordList;			
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,3));
			if ~isempty(sppsEmptyMinorValence)   
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve(seedPoints_(...
						sppsEmptyMinorValence,:), minorPSL.phyCoordList);					
				potentialSolidSppsMinor = find(potentialDisListMinor<relaxedFactor_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,3) = 1;				
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');					
					seedPointsValence_(modifiedMajorValences,1) = 1;
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValence_(modifiedMediumValences,2) = 1;					
				end
			end					
		end
		looper = sum(sum(seedPointsValence_));
		disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
			' Total.: ' sprintf('%6i',3*numSeedPoints)]);			
	end
	PSLs_CompactPSLs();
end

function [potentialDisList, potentialPosList] = GetDisListOfPointList2Curve(pointList, curveLine)
	global mergeTrigger_;	
	disT = (curveLine(:,end-2) - pointList(:,end-2)').^2;
	disT = disT + (curveLine(:,end-1) - pointList(:,end-1)').^2;
	disT = disT + (curveLine(:,end) - pointList(:,end)').^2;
	disT = sqrt(disT);	
	[minVal, minValPos] = min(disT,[],1);
	potentialDisList = minVal';
	potentialDisList = potentialDisList/mergeTrigger_;
	potentialPosList = curveLine(minValPos,:);
end

function modifiedValences = HighCurvatureModification(spps2BeMerged, psDir)
	global majorCoordList_; 
	global mediumCoordList_; 
	global minorCoordList_;		
	global seedPoints_;
	global seedPointsValence_;
	global mergeTrigger_;
	global relaxedFactor_;

	coordList = [];
	switch psDir
		case 'MAJOR'
			if isempty(majorCoordList_), modifiedValences = []; return; end
			coordList = majorCoordList_;
			% spps2BeMerged = spps2BeMerged(find(0==seedPointsValence_(spps2BeMerged,1)));
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,1));
		case 'MEDIUM'
			if isempty(mediumCoordList_), modifiedValences = []; return; end
			coordList = mediumCoordList_;
			% spps2BeMerged = spps2BeMerged(find(0==seedPointsValence_(spps2BeMerged,2)));
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,2));
		case 'MINOR'
			if isempty(minorCoordList_), modifiedValences = []; return; end
			coordList = minorCoordList_;
			% spps2BeMerged = spps2BeMerged(find(0==seedPointsValence_(spps2BeMerged,3)));
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,3));
	end
	pointList = seedPoints_(spps2BeMerged,end-2:end);
	disT = (coordList(:,1) - pointList(:,1)').^2;
	disT = disT + (coordList(:,2) - pointList(:,2)').^2;
	disT = disT + (coordList(:,3) - pointList(:,3)').^2;
	disT = sqrt(disT);		
	minVal = min(disT, [], 1);
	minVal = minVal/mergeTrigger_;
	switch psDir
		case 'MAJOR'
			modifiedValences = find(minVal<relaxedFactor_);	
		case 'MEDIUM'
			modifiedValences = find(minVal<relaxedFactor_);	
		case 'MINOR'
			modifiedValences = find(minVal<relaxedFactor_);	
	end	
	modifiedValences = spps2BeMerged(modifiedValences);
end

function PSLs_CompactPSLs()
	global minLengthVisiblePSLs_;
	global majorPSLpool_; 
	global mediumPSLpool_; 
	global minorPSLpool_;
	
	filterThreshold = minLengthVisiblePSLs_;
	numMajorPSLs = length(majorPSLpool_);
	tarIndice = [];
	for ii=1:numMajorPSLs
		if majorPSLpool_(ii).length > filterThreshold
			tarIndice(end+1,1) = ii;
		end
	end
	majorPSLpool_ = majorPSLpool_(tarIndice);

	numMediumPSLs = length(mediumPSLpool_);
	tarIndice = [];
	for ii=1:numMediumPSLs
		if mediumPSLpool_(ii).length > filterThreshold
			tarIndice(end+1,1) = ii;
		end
	end
	mediumPSLpool_ = mediumPSLpool_(tarIndice);	

	numMinorPSLs = length(minorPSLpool_);
	tarIndice = [];
	for ii=1:numMinorPSLs
		if minorPSLpool_(ii).length > filterThreshold
			tarIndice(end+1,1) = ii;
		end
	end
	minorPSLpool_ = minorPSLpool_(tarIndice);
end

function iPSL = CreatePrincipalStressLine(startPoint, tracingType)
	global integrationStepLimit_;
	iPSL = Data_PrincipalStressLineStruct();
	switch tracingType
		case 'MAJOR', psDir = [10 11 12];
		case 'MEDIUM', psDir = [6 7 8];		
		case 'MINOR', psDir = [2 3 4];
	end
	%%1. prepare for tracing			
	[eleIndex, cartesianStress, principalStress, opt] = PSLs_PreparingForTracing(startPoint);
	if 0==opt, return; end
	
	%%2. tracing PSL
	startPoint = startPoint(end-2:end);
	PSLphyCoordList = startPoint;
	PSLcartesianStressList = cartesianStress;
	PSLeleIndexList = eleIndex;
	PSLprincipalStressList = principalStress;			
	%%2.1 along first direction (v1)		
	[phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
		TracingPSL_RK2_CartesianMesh(startPoint, principalStress(1,psDir), eleIndex, psDir, integrationStepLimit_);		
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	PSLcartesianStressList = [PSLcartesianStressList; cartesianStressList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLprincipalStressList = [PSLprincipalStressList; principalStressList];
	%%2.2 along second direction (-v1)	
	[phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
		TracingPSL_RK2_CartesianMesh(startPoint, -principalStress(1,psDir), eleIndex, psDir, integrationStepLimit_);		
	if size(phyCoordList,1) > 1
		phyCoordList = flip(phyCoordList);
		cartesianStressList = flip(cartesianStressList);
		eleIndexList = flip(eleIndexList);
		principalStressList = flip(principalStressList);
	end						
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	PSLcartesianStressList = [cartesianStressList; PSLcartesianStressList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLprincipalStressList = [principalStressList; PSLprincipalStressList];
	%%2.3 finish Tracing the current major PSL	
	iPSL.midPointPosition = size(phyCoordList,1)+1;
	iPSL.length = size(PSLphyCoordList,1);
	iPSL.eleIndexList = PSLeleIndexList;
	iPSL.phyCoordList = PSLphyCoordList;
	iPSL.cartesianStressList = PSLcartesianStressList;	
	iPSL.principalStressList = PSLprincipalStressList;	
end

function [eleIndex, cartesianStress, principalStress, opt] = PSLs_PreparingForTracing(startPoint)
	global nodeCoords_; 
	global meshHierarchy_; 
	global nodStruct_;
	global fieldDataType_;
	global cartesianStressField_;
	global frameField_;
	global eleCentroidList_;
	global meshType_;

	eleIndex = 0;
	cartesianStress = 0;
	principalStress = 0;	
	[targetEleIndex, paraCoordinates, opt] = PSLs_SearchNextIntegratingPointOnCartesianMesh(startPoint(end-2:end));	
	if 0==opt, return; end
	eleIndex = double(targetEleIndex);
	NIdx = meshHierarchy_(1).eNodMat(eleIndex,:)';
	eleCartesianStress = cartesianStressField_(NIdx,:);				
	cartesianStress = ...
		FEA_ShapeFunction(paraCoordinates(1), paraCoordinates(2), paraCoordinates(3)) * eleCartesianStress;		
	principalStress = FEA_ComputePrincipalStress(cartesianStress);	
end

function [phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
			TracingPSL_RK2_CartesianMesh(startPoint, iniDir, elementIndex, typePSL, limiSteps)
	global meshHierarchy_;
	global nodeCoords_;
	global cartesianStressField_;
	global tracingStepWidth_;

	phyCoordList = zeros(limiSteps,3);
	cartesianStressList = zeros(limiSteps,6);
	eleIndexList = zeros(limiSteps,1);
	principalStressList = zeros(limiSteps,12);	

	%%initialize initial k1 and k2
	k1 = iniDir;
	midPot = startPoint + k1*tracingStepWidth_/2;
	index = 0;	
	[elementIndex, paraCoordinates, bool1] = PSLs_SearchNextIntegratingPointOnCartesianMesh(midPot);		
	if bool1
		cartesianStress = cartesianStressField_(meshHierarchy_(1).eNodMat(elementIndex,:)', :);
		cartesianStressOnGivenPoint = ...
			FEA_ShapeFunction(paraCoordinates(1), paraCoordinates(2), paraCoordinates(3)) * cartesianStress;		
		principalStress = FEA_ComputePrincipalStress(cartesianStressOnGivenPoint);		
 		[k2, terminationCond] = PSLs_BidirectionalFeatureProcessing(k1, principalStress(typePSL));
		nextPoint = startPoint + tracingStepWidth_*k2;
		[elementIndex, paraCoordinates, bool1] = PSLs_SearchNextIntegratingPointOnCartesianMesh(nextPoint);
		while bool1
			index = index + 1; if index > limiSteps, index = index-1; break; end
			%%k1
			cartesianStress = cartesianStressField_(meshHierarchy_(1).eNodMat(elementIndex,:)', :);
			cartesianStressOnGivenPoint = ...
				FEA_ShapeFunction(paraCoordinates(1), paraCoordinates(2), paraCoordinates(3)) * cartesianStress;
			principalStress = FEA_ComputePrincipalStress(cartesianStressOnGivenPoint);
			evs = principalStress([1 5 9]);		
 			[k1, terminationCond] = PSLs_BidirectionalFeatureProcessing(iniDir, principalStress(typePSL));		
            if ~terminationCond, index = index-1; break; end
			%%k2
			midPot = nextPoint + k1*tracingStepWidth_/2;
			[elementIndex2, paraCoordinates2, bool1] = PSLs_SearchNextIntegratingPointOnCartesianMesh(midPot);
			if ~bool1, index = index-1; break; end
			cartesianStress2 = cartesianStressField_(meshHierarchy_(1).eNodMat(elementIndex2,:)', :);
			cartesianStressOnGivenPoint2 = ...
				FEA_ShapeFunction(paraCoordinates2(1), paraCoordinates2(2), paraCoordinates2(3)) * cartesianStress2;
			principalStress2 = FEA_ComputePrincipalStress(cartesianStressOnGivenPoint2);
 			[k2, terminationCond] = PSLs_BidirectionalFeatureProcessing(k1, principalStress2(typePSL));
			%%store	
			iniDir = k1;
			phyCoordList(index,:) = nextPoint;
			cartesianStressList(index,:) = cartesianStressOnGivenPoint;
			eleIndexList(index,:) = elementIndex;			
			principalStressList(index,:) = principalStress;
			%%next point
			nextPoint = nextPoint + tracingStepWidth_*k2;
			[elementIndex, paraCoordinates, bool1] = PSLs_SearchNextIntegratingPointOnCartesianMesh(nextPoint);				
		end
	end
	phyCoordList = phyCoordList(1:index,:);
	cartesianStressList = cartesianStressList(1:index,:);
	eleIndexList = eleIndexList(1:index,:);
	principalStressList = principalStressList(1:index,:);	
end

function val = Data_PrincipalStressLineStruct()
	val = struct(...
		'length',						0,	...
		'midPointPosition',				0,	...		
		'phyCoordList',					[], ...
		'eleIndexList',					[], ...
		'adjacentVoxels',				[],	...
		'paraCoordList',				[], ...
		'cartesianStressList',			[],	...
		'vonMisesStressList',			[], ...
		'principalStressList', 			[], ...
		'tubePatchIndices',				[],	...
		'strainEnergyApprox',			0,	...
		'RibbonPatchIndices',			[],	...
		'RibbonOutlinePatchIndices',	[]	...		
	);	
end

function [targetDirection, terminationCond] = PSLs_BidirectionalFeatureProcessing(originalVec, Vec)
	global permittedMaxAdjacentTangentAngleDeviation_;
	terminationCond = 1;
	angle1 = acos(originalVec*Vec');
	angle2 = acos(-originalVec*Vec');	
	if angle1 < angle2
		targetDirection = Vec;
		if angle1 > pi/permittedMaxAdjacentTangentAngleDeviation_, terminationCond = 0; end
	else
		targetDirection = -Vec;
		if angle2 > pi/permittedMaxAdjacentTangentAngleDeviation_, terminationCond = 0; end
	end
end

function [nextElementIndex, paraCoordinates, opt] = PSLs_SearchNextIntegratingPointOnCartesianMesh(physicalCoordinates)
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
		relatedNodes = meshHierarchy_(1).eNodMat(nextElementIndex,:);
		relatedNodeCoords = nodeCoords_(relatedNodes',:)-boundingBox_(1,:);
		paraCoordinates = 2*(physicalCoordinates - relatedNodeCoords(1,:)) / meshHierarchy_(1).eleSize(1) - 1;
	end	
end