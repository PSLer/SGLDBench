%% This code is created for performing the 3D-TSV to the stress tensor field simulated on the 1st-order Quadrilateral and Triangular Mesh 
% Easy-Use Version
% Author: Junpeng Wang (junpeng.wang@tum.de)
% Date: 2023-12-12
% clear all; clc;

% global majorPSLpoolSmoothed_;
% global mediumPSLpoolSmoothed_; 
% global minorPSLpoolSmoothed_; 


%%1. Import Data
% stressfileName = 'StressField_Tet_v2_eleWise.stress';
% ImportStressFields(stressfileName);


%%2. PSLs generation
% PSLsDensityCtrl = 10; %%The larger this value is, the denser the PSLs are.
% TSV3D(PSLsDensityCtrl,10);
% figure; ShowPSLs(100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TSV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EZ_TSV3D_eleWise(PSLsDensityCtrl, varargin)
	global outPath_;
	global boundingBoxSmoothed_;
	global nodeCoordsSmoothed_;
	global eleCentroidListSmoothed_;
	global numElesSmoothed_;
	global boundaryElementsSmoothed_;
	global eleSizeListSmoothed_;
	
	global mergingThresholdSmoothed_;
	global mergingThresholdMapSmoothed_;
	global tracingStepWidthSmoothed_;
	global integrationStepLimitSmoothed_;
	global permittedMaxAdjacentTangentAngleDeviationSmoothed_;
	global relaxedFactorSmoothed_;
	
	global seedPointsHistorySmoothed_;
	global seedAssociatedElesSmoothed_;
	global seedPointsSmoothed_;
	global seedPointsValenceSmoothed_;	
	global majorPSLpoolSmoothed_;
	global mediumPSLpoolSmoothed_;
    global minorPSLpoolSmoothed_; 
	global majorCoordListSmoothed_;
	global mediumCoordListSmoothed_; 
    global minorCoordListSmoothed_;

	%%Import Data
	stressfileName = strcat(outPath_, 'StressField_Tet_v2_eleWise.stress');
	ImportStressFields(stressfileName);	
	
	%%Settings
	mergingThresholdSmoothed_ = min(boundingBoxSmoothed_(2,:)-boundingBoxSmoothed_(1,:))/PSLsDensityCtrl;
	permittedMaxAdjacentTangentAngleDeviationSmoothed_ = 10;
	tracingStepWidthSmoothed_ = 1.0 * min(boundingBoxSmoothed_(2,:)-boundingBoxSmoothed_(1,:))/100;
	tracingStepWidthSmoothed_ = 0.5 * eleSizeListSmoothed_;
	integrationStepLimitSmoothed_ = ceil(1.5*norm(boundingBoxSmoothed_(2,:)-boundingBoxSmoothed_(1,:))/median(tracingStepWidthSmoothed_));
	relaxedFactorSmoothed_ = 1.0;
	
    if 1==nargin
        seedDensityCtrl = 1;
    else
        seedDensityCtrl = varargin{1};
    end
	seedAssociatedElesSmoothed_ = 1:seedDensityCtrl:numElesSmoothed_; seedAssociatedElesSmoothed_ = seedAssociatedElesSmoothed_(:);
	% seedAssociatedElesSmoothed_(boundaryElementsSmoothed_) = [];
	seedPointsHistorySmoothed_ = eleCentroidListSmoothed_(seedAssociatedElesSmoothed_,:);
	seedPointsValenceSmoothed_ = ones(size(seedPointsHistorySmoothed_));
	numSeedPoints = size(seedPointsHistorySmoothed_,1);
	seedPointsSmoothed_ = [seedAssociatedElesSmoothed_ seedPointsHistorySmoothed_];
	
	InitializeMergingThresholdMap();
	
	majorPSLpoolSmoothed_ = PrincipalStressLineStruct();
	mediumPSLpoolSmoothed_ = PrincipalStressLineStruct();
	minorPSLpoolSmoothed_ = PrincipalStressLineStruct();

	%% Exclude the irrelated Principal Stress Fields
	selectedPrincipalStressField_ = [1 2 3];
	numPSF = length(selectedPrincipalStressField_);
	for ii=1:numPSF
		iPSF = selectedPrincipalStressField_(ii);
		switch iPSF
			case 1, seedPointsValenceSmoothed_(:,1) = 0; 
			case 2, seedPointsValenceSmoothed_(:,2) = 0; 
			case 3, seedPointsValenceSmoothed_(:,3) = 0;
		end
	end

	PreprocessSeedPoints();
	
	majorCoordListSmoothed_ = [];
	mediumCoordListSmoothed_ = [];
	minorCoordListSmoothed_ = [];
	
	startCoord_ = boundingBoxSmoothed_(1,:) + sum(boundingBoxSmoothed_, 1)/2;
	
	%% Seeding
	its = 0;
	looper = sum(sum(seedPointsValenceSmoothed_));
	while looper<3*numSeedPoints
		its = its + 1;
		valenceMetric = sum(seedPointsValenceSmoothed_,2);
		%% 1st Priority: semi-empty seeds > empty seeds, which helps get PSLs intersection
		%% 2nd Priority: seeds with same valence, the one closest to the start point goes first	
		switch numPSF
			case 1
				unFinishedSppsValence2 = find(2==valenceMetric);
				[~, tarPos] = min(vecnorm(startCoord_-seedPointsSmoothed_(unFinishedSppsValence2,end-2:end),2,2));
				spp = unFinishedSppsValence2(tarPos);
			case 2
				unFinishedSppsValence2 = find(2==valenceMetric); 
				if ~isempty(unFinishedSppsValence2) %% 1st Priority
					[~, tarPos] = min(vecnorm(startCoord_-seedPointsSmoothed_(unFinishedSppsValence2,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence2(tarPos);
				else
					unFinishedSppsValence1 = find(1==valenceMetric);
					[~, tarPos] = min(vecnorm(startCoord_-seedPointsSmoothed_(unFinishedSppsValence1,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence1(tarPos);			
				end					
			case 3
				unFinishedSppsValence12 = find(3>valenceMetric); 
				unFinishedSppsValence12 = unFinishedSppsValence12(valenceMetric(unFinishedSppsValence12)>0); 
				if ~isempty(unFinishedSppsValence12) %% 1st Priority
					[~, tarPos] = min(vecnorm(startCoord_-seedPointsSmoothed_(unFinishedSppsValence12,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence12(tarPos);
				else
					unFinishedSppsValence0 = find(0==valenceMetric);
					[~, tarPos] = min(vecnorm(startCoord_-seedPointsSmoothed_(unFinishedSppsValence0,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence0(tarPos);		
				end						
		end
		valences = seedPointsValenceSmoothed_(spp,:);						
		seed = seedPointsSmoothed_(spp,:);		
		
		if 0==valences(1)
			seedPointsValenceSmoothed_(spp,1) = 1;
			majorPSL = Have1morePSL(seed, 'MAJOR');		
			if 0==majorPSL.length
				looper = sum(sum(seedPointsValenceSmoothed_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end			
			majorPSLpoolSmoothed_(end+1,1) = majorPSL;				
			majorCoordListSmoothed_(end+1:end+majorPSL.length,:) = majorPSL.phyCoordList;
			sppsEmptyMajorValence = find(0==seedPointsValenceSmoothed_(:,1));
			if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve(seedPointsSmoothed_(...
						sppsEmptyMajorValence,:), [majorPSL.eleIndexList majorPSL.phyCoordList], 'MAJOR');					
				potentialSolidSppsMajor = find(potentialDisListMajor<relaxedFactorSmoothed_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);
					seedPointsSmoothed_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);								
					seedPointsValenceSmoothed_(spps2BeMerged,1) = 1;
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValenceSmoothed_(modifiedMediumValences,2) = 1;					
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');				
					seedPointsValenceSmoothed_(modifiedMinorValences,3) = 1;	
				end
			end				
		end		

		if 0==valences(2)
			seedPointsValenceSmoothed_(spp,2) = 1;
			mediumPSL = Have1morePSL(seed, 'MEDIUM');		
			if 0==mediumPSL.length
				looper = sum(sum(seedPointsValenceSmoothed_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end
			mediumPSLpoolSmoothed_(end+1,1) = mediumPSL;
			mediumCoordListSmoothed_(end+1:end+mediumPSL.length,:) = mediumPSL.phyCoordList;
			sppsEmptyMediumValence = find(0==seedPointsValenceSmoothed_(:,2));
			if ~isempty(sppsEmptyMediumValence)
				[potentialDisListMedium, potentialPosListMedium] = GetDisListOfPointList2Curve(seedPointsSmoothed_(...
						sppsEmptyMediumValence,:), [mediumPSL.eleIndexList mediumPSL.phyCoordList], 'MEDIUM');					
				potentialSolidSppsMedium = find(potentialDisListMedium<relaxedFactorSmoothed_);
				if ~isempty(potentialSolidSppsMedium)
					spps2BeMerged = sppsEmptyMediumValence(potentialSolidSppsMedium);
					seedPointsSmoothed_(spps2BeMerged,:) = potentialPosListMedium(potentialSolidSppsMedium,:);								
					seedPointsValenceSmoothed_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValenceSmoothed_(modifiedMajorValences,1) = 1;
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValenceSmoothed_(modifiedMinorValences,3) = 1;
				end
			end				
		end

		if 0==valences(3)
			seedPointsValenceSmoothed_(spp,3) = 1;			
			minorPSL = Have1morePSL(seed, 'MINOR');			
			if 0==minorPSL.length
				looper = sum(sum(seedPointsValenceSmoothed_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end		
			minorPSLpoolSmoothed_(end+1,1) = minorPSL;
			minorCoordListSmoothed_(end+1:end+minorPSL.length,:) = minorPSL.phyCoordList;		
			sppsEmptyMinorValence = find(0==seedPointsValenceSmoothed_(:,3));
			if ~isempty(sppsEmptyMinorValence)   
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve(seedPointsSmoothed_(...
						sppsEmptyMinorValence,:), [minorPSL.eleIndexList minorPSL.phyCoordList], 'MINOR');					
				potentialSolidSppsMinor = find(potentialDisListMinor<relaxedFactorSmoothed_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPointsSmoothed_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValenceSmoothed_(spps2BeMerged,3) = 1;				
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');					
					seedPointsValenceSmoothed_(modifiedMajorValences,1) = 1;
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValenceSmoothed_(modifiedMediumValences,2) = 1;					
				end
			end					
		end		
		looper = sum(sum(seedPointsValenceSmoothed_));
		disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
			' Total.: ' sprintf('%6i',3*numSeedPoints)]);			
	end

	majorPSLpoolSmoothed_ = CompactStreamlines(majorPSLpoolSmoothed_, 5);
	mediumPSLpoolSmoothed_ = CompactStreamlines(mediumPSLpoolSmoothed_, 5);
	minorPSLpoolSmoothed_ = CompactStreamlines(minorPSLpoolSmoothed_, 5);	
end

function InitializeMergingThresholdMap()
	global numElesSmoothed_;
	global mergingThresholdSmoothed_;
	global mergingThresholdMapSmoothed_;
	mergingThresholdMapSmoothed_ = repmat(mergingThresholdSmoothed_, numElesSmoothed_, 3);
end

function PreprocessSeedPoints()
	global seedPointsSmoothed_;
	global seedPointsValenceSmoothed_;
	global majorPSLpoolSmoothed_; 
	global mediumPSLpoolSmoothed_; 
	global minorPSLpoolSmoothed_; 
	global relaxedFactorSmoothed_;
	global multiMergingThresholdsCtrlSmoothed_;
	
	numMajorPSLs = length(majorPSLpoolSmoothed_);
	for ii=1:numMajorPSLs
		majorPSL = majorPSLpoolSmoothed_(ii);
		if majorPSL.length>0					
			sppsEmptyMajorValence = find(0==seedPointsValenceSmoothed_(:,1));
            if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve(...	
					seedPointsSmoothed_(sppsEmptyMajorValence,:), [majorPSL.eleIndexList majorPSL.phyCoordList], 'MAJOR');
				potentialSolidSppsMajor = find(potentialDisListMajor<=multiMergingThresholdsCtrlSmoothed_(1)*relaxedFactorSmoothed_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);							
					seedPointsSmoothed_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);
					seedPointsValenceSmoothed_(spps2BeMerged,1) = 1;						
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValenceSmoothed_(modifiedMediumValences,2) = 1;					
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValenceSmoothed_(modifiedMinorValences,3) = 1;							
				end
			end
		end
	end

	numMediumPSLs = length(mediumPSLpoolSmoothed_);
	for ii=1:numMediumPSLs
		mediumPSL = mediumPSLpoolSmoothed_(ii);
		if mediumPSL.length>0					
			sppsEmptyMediumValence = find(0==seedPointsValenceSmoothed_(:,2));
            if ~isempty(sppsEmptyMediumValence)
				[potentialDisListMedium, potentialPosListMedium] = GetDisListOfPointList2Curve(...	
					seedPointsSmoothed_(sppsEmptyMediumValence,:), [mediumPSL.eleIndexList mediumPSL.phyCoordList], 'MEDIUM');
				potentialSolidSppsMedium = find(potentialDisListMedium<=multiMergingThresholdsCtrlSmoothed_(2)*relaxedFactorSmoothed_);
				if ~isempty(potentialSolidSppsMedium)
					spps2BeMerged = sppsEmptyMediumValence(potentialSolidSppsMedium);							
					seedPointsSmoothed_(spps2BeMerged,:) = potentialPosListMedium(potentialSolidSppsMedium,:);
					seedPointsValenceSmoothed_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValenceSmoothed_(modifiedMajorValences,1) = 1;						
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValenceSmoothed_(modifiedMinorValences,3) = 1;							
				end
			end
		end
	end

	numMinorPSLs = length(minorPSLpoolSmoothed_);
	for ii=1:numMinorPSLs
		minorPSL = minorPSLpoolSmoothed_(ii);
		if minorPSL.length>0	
			sppsEmptyMinorValence = find(0==seedPointsValenceSmoothed_(:,3));
            if ~isempty(sppsEmptyMinorValence)
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve(...	
					seedPointsSmoothed_(sppsEmptyMinorValence,:), [minorPSL.eleIndexList minorPSL.phyCoordList], 'MINOR');
				potentialSolidSppsMinor = find(potentialDisListMinor<=multiMergingThresholdsCtrlSmoothed_(3)*relaxedFactorSmoothed_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPointsSmoothed_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValenceSmoothed_(spps2BeMerged,3) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValenceSmoothed_(modifiedMajorValences,1) = 1;	
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValenceSmoothed_(modifiedMediumValences,2) = 1;						
				end
			end
		end
	end		
end

function modifiedValences = HighCurvatureModification(spps2BeMerged, psDir)
	global majorCoordListSmoothed_;
    global mediumCoordListSmoothed_;
	global minorCoordListSmoothed_;		
	global seedPointsSmoothed_;
	global seedPointsValenceSmoothed_;
	global mergingThresholdSmoothed_;
	global relaxedFactorSmoothed_;

	coordList = [];
	switch psDir
		case 'MAJOR'
			if isempty(majorCoordListSmoothed_), modifiedValences = []; return; end
			coordList = majorCoordListSmoothed_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValenceSmoothed_(spps2BeMerged,1));
		case 'MEDIUM'
			if isempty(mediumCoordListSmoothed_), modifiedValences = []; return; end
			coordList = mediumCoordListSmoothed_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValenceSmoothed_(spps2BeMerged,2));			
		case 'MINOR'
			if isempty(minorCoordListSmoothed_), modifiedValences = []; return; end
			coordList = minorCoordListSmoothed_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValenceSmoothed_(spps2BeMerged,3));
	end
	pointList = seedPointsSmoothed_(spps2BeMerged,end-2:end);
	disT = (coordList(:,1) - pointList(:,1)').^2;
	disT = disT + (coordList(:,2) - pointList(:,2)').^2;
	disT = disT + (coordList(:,3) - pointList(:,3)').^2;
	disT = sqrt(disT);		
	minVal = min(disT, [], 1);
	minVal = minVal/mergingThresholdSmoothed_;
	modifiedValences = find(minVal<relaxedFactorSmoothed_);	
	modifiedValences = spps2BeMerged(modifiedValences);
end

function [potentialDisList, potentialPosList] = GetDisListOfPointList2Curve(pointList, curveLine, psDir)
	global mergingThresholdMapSmoothed_;
	global mergingThresholdSmoothed_;	
	disT = (curveLine(:,end-2) - pointList(:,end-2)').^2;
	disT = disT + (curveLine(:,end-1) - pointList(:,end-1)').^2;
	disT = disT + (curveLine(:,end) - pointList(:,end)').^2;
	disT = sqrt(disT);	
	[minVal, minValPos] = min(disT,[],1);

	potentialPosList = curveLine(minValPos,:);
	switch psDir
		case 'MAJOR', idx = 1;
		case 'MEDIUM', idx = 2;
		case 'MINOR', idx = 3;
	end
	
	potentialDisList = minVal';
	% potentialDisList = potentialDisList/mergingThresholdSmoothed_;
	potentialDisList = potentialDisList ./ mergingThresholdMapSmoothed_(potentialPosList(:,1),idx);
end

function [eleIndex, cartesianStress, principalStress, opt] = PreparingForTracing(startPoint)
	global nodeCoordsSmoothed_; 
	global eNodMatSmoothed_; 
	global nodStructSmoothed_;
	global fieldDataTypeSmoothed_;
	global principalStressFieldSmoothed_;
	global frameFieldSmoothed_;
	global eleCentroidListSmoothed_;
	global meshTypeSmoothed_;

	eleIndex = 0;
	cartesianStress = 0;
	principalStress = 0;	
	switch numel(startPoint)
		case 3
			disList = vecnorm(startPoint-eleCentroidListSmoothed_, 2, 2);
			[~, targetEleIndex0] = min(disList);	
			[eleIndex, opt] = PositioningOnUnstructuredMesh(targetEleIndex0, startPoint);
			if ~opt, return; end			
		case 4
			eleIndex = startPoint(1);
			startPoint = startPoint(2:4);
            opt = 1;			
		otherwise
			error('Wrong Input For the Seed!')
	end
	NIdx = eNodMatSmoothed_(eleIndex,:)';
	eleNodeCoords = nodeCoordsSmoothed_(NIdx,:);
	principalStress = principalStressFieldSmoothed_(eleIndex,:);
end

function iPSL = Have1morePSL(startPoint, tracingType)
	global integrationStepLimitSmoothed_;

	iPSL = PrincipalStressLineStruct();
	switch tracingType
		case 'MAJOR', psDir = [10 11 12];
		case 'MEDIUM', psDir = [6 7 8];
		case 'MINOR', psDir = [2 3 4];
	end

	%%1. prepare for tracing			
	[eleIndex, cartesianStress, principalStress, opt] = PreparingForTracing(startPoint);
	if 0==opt, return; end
	
	%%2. tracing PSL
	startPoint = startPoint(end-2:end);
	PSLphyCoordList = startPoint;
	% PSLcartesianStressList = cartesianStress;
	PSLeleIndexList = eleIndex;
	PSLprincipalStressList = principalStress;
	
	%%2.1 along first direction (v1)		
	[phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
		TracingPSL_RK2(startPoint, principalStress(1,psDir), eleIndex, psDir, integrationStepLimitSmoothed_);		
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	% PSLcartesianStressList = [PSLcartesianStressList; cartesianStressList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLprincipalStressList = [PSLprincipalStressList; principalStressList];
	
	%%2.2 along second direction (-v1)	
	[phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
		TracingPSL_RK2(startPoint, -principalStress(1,psDir), eleIndex, psDir, integrationStepLimitSmoothed_);		
	if size(phyCoordList,1) > 1
		phyCoordList = flip(phyCoordList);
		cartesianStressList = flip(cartesianStressList);
		eleIndexList = flip(eleIndexList);
		principalStressList = flip(principalStressList);
	end						
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	% PSLcartesianStressList = [cartesianStressList; PSLcartesianStressList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLprincipalStressList = [principalStressList; PSLprincipalStressList];
	
	%%2.3 finish Tracing the current major PSL	
	iPSL.midPointPosition = size(phyCoordList,1)+1;
	iPSL.length = size(PSLphyCoordList,1);
	iPSL.eleIndexList = PSLeleIndexList;
	iPSL.phyCoordList = PSLphyCoordList;
	% iPSL.cartesianStressList = PSLcartesianStressList;	
	iPSL.principalStressList = PSLprincipalStressList;		
end

function val = PrincipalStressLineStruct()
	val = struct(...
		'ith',						0, 	...
		'length',					0,	...
		'midPointPosition',			0,	...		
		'phyCoordList',				[], ...
		'eleIndexList',				[], ...
		'cartesianStressList',		[],	...
		'principalStressList',		[] ...
	);	
end

function [phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
			TracingPSL_RK2(startPoint, iniDir, elementIndex, typePSL, limiSteps)			
	global eNodMatSmoothed_;
	global nodeCoordsSmoothed_;
	% global cartesianStressFieldSmoothed_;
	global principalStressFieldSmoothed_;
	global tracingStepWidthSmoothed_;
	phyCoordList = zeros(limiSteps,3);
	cartesianStressList = zeros(limiSteps,6);
	eleIndexList = zeros(limiSteps,1);
	principalStressList = zeros(limiSteps,12);	
	
	index = 0;
	k1 = iniDir;
	%%re-scale stepsize if necessary
	stepsize = tracingStepWidthSmoothed_(elementIndex);
	testPot = startPoint + k1*tracingStepWidthSmoothed_(elementIndex);
	[~, testPot1, bool1] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, testPot, startPoint, 1);
	if bool1
		stepsize = norm(testPot1-startPoint)/norm(testPot-startPoint) * stepsize;
		%%initialize initial k1 and k2
		midPot = startPoint + k1*stepsize/2;
		[elementIndex2, ~, bool2] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, midPot, startPoint, 0);
		if bool2 %%just in case
			% NIdx = eNodMatSmoothed_(elementIndex2,:)';
			% vtxStress = cartesianStressFieldSmoothed_(NIdx, :);
			% vtxCoords = nodeCoordsSmoothed_(NIdx,:);
			% cartesianStressOnGivenPoint = ElementInterpolationInverseDistanceWeighting(vtxCoords, vtxStress, midPot);
			% principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);
			principalStress = principalStressFieldSmoothed_(elementIndex2,:);
			[k2, ~] = BidirectionalFeatureProcessing(k1, principalStress(typePSL));
			nextPoint = startPoint + stepsize*k2;
			[elementIndex, ~, bool3] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, nextPoint, startPoint, 0);
			while bool3
				index = index + 1; if index > limiSteps, index = index-1; break; end
				% NIdx = eNodMatSmoothed_(elementIndex,:)';
				% vtxStress = cartesianStressFieldSmoothed_(NIdx, :);
				% vtxCoords = nodeCoordsSmoothed_(NIdx,:); 
				% cartesianStressOnGivenPoint = ElementInterpolationInverseDistanceWeighting(vtxCoords, vtxStress, nextPoint); 
				% principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);
				principalStress = principalStressFieldSmoothed_(elementIndex,:);				
				%%k1
				[k1, terminationCond] = BidirectionalFeatureProcessing(iniDir, principalStress(typePSL));	
				if ~terminationCond, index = index-1; break; end
				%%k2
				%%re-scale stepsize if necessary
				stepsize = tracingStepWidthSmoothed_(elementIndex);
				testPot = nextPoint + k1*stepsize;
				[~, testPot1, bool1] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, testPot, nextPoint, 1);			
				if ~bool1, index = index-1; break; end
				stepsize = norm(testPot1-nextPoint)/norm(testPot-nextPoint) * stepsize;
				midPot = nextPoint + k1*stepsize/2;
				[elementIndex2, ~, bool1] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, midPot, nextPoint, 0);					
				if ~bool1, index = index-1; break; end	
				% NIdx2 = eNodMatSmoothed_(elementIndex2,:)';	
				% vtxStress2 = cartesianStressFieldSmoothed_(NIdx2,:);
				% vtxCoords2 = nodeCoordsSmoothed_(NIdx2,:);
				% cartesianStressOnGivenPoint2 = ElementInterpolationInverseDistanceWeighting(vtxCoords2, vtxStress2, midPot);
				% principalStress2 = ComputePrincipalStress(cartesianStressOnGivenPoint2);
				principalStress2 = principalStressFieldSmoothed_(elementIndex2,:);	
				[k2, ~] = BidirectionalFeatureProcessing(k1, principalStress2(typePSL));					
				%%store	
				iniDir = k1;
				phyCoordList(index,:) = nextPoint;
				% cartesianStressList(index,:) = cartesianStressOnGivenPoint;
				eleIndexList(index,:) = elementIndex;
				principalStressList(index,:) = principalStress;
				%%next point
				nextPoint0 = nextPoint + stepsize*k2;
				[elementIndex, ~, bool3] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, nextPoint0, nextPoint, 0);			
				nextPoint = nextPoint0;				
			end
		end
	end
	phyCoordList = phyCoordList(1:index,:);
	% cartesianStressList = cartesianStressList(1:index,:);
	eleIndexList = eleIndexList(1:index,:);
	principalStressList = principalStressList(1:index,:);	
end

function [targetDirection, terminationCond] = BidirectionalFeatureProcessing(originalVec, Vec)
	global permittedMaxAdjacentTangentAngleDeviationSmoothed_;
	terminationCond = 1;
	angle1 = acos(originalVec*Vec');
	angle2 = acos(-originalVec*Vec');	
	if angle1 < angle2
		targetDirection = Vec;
		if angle1 > pi/permittedMaxAdjacentTangentAngleDeviationSmoothed_, terminationCond = 0; end
	else
		targetDirection = -Vec;
		if angle2 > pi/permittedMaxAdjacentTangentAngleDeviationSmoothed_, terminationCond = 0; end
	end	
end

function [eleIndex, opt] = PositioningOnUnstructuredMesh(targetEleIndex0, startPoint)
	global eNodMatSmoothed_; 
	global nodStructSmoothed_;
	global eleCentroidListSmoothed_;
	opt = IsThisPointWithinThatElement(targetEleIndex0, startPoint);
	if opt
		eleIndex = targetEleIndex0;		
	else %% Search the Adjacent Elements
		tarNodes = eNodMatSmoothed_(targetEleIndex0,:);
		allPotentialAdjacentElements = unique([nodStructSmoothed_(tarNodes(:)).adjacentEles]);
		potentialAdjacentElements = setdiff(allPotentialAdjacentElements, targetEleIndex0);
		for ii=1:length(potentialAdjacentElements)
			iEle = potentialAdjacentElements(ii);
			opt = IsThisPointWithinThatElement(iEle, startPoint);
			if opt, eleIndex = iEle; break; end
		end
	end
	if 0==opt		
		disList = vecnorm(startPoint-eleCentroidListSmoothed_(allPotentialAdjacentElements,:), 2, 2);
		[~, nearOptimalEle] = min(disList);
		eleIndex = allPotentialAdjacentElements(nearOptimalEle);
	end
end

function opt = IsThisPointWithinThatElement(tarEleIndex, iCoord)
	global eleStructSmoothed_;
	global nodeCoordsSmoothed_; 
	global eNodMatSmoothed_;
	global numFacesPerEleSmoothed_;
	opt = 1; 
	iElefaceCentres = eleStructSmoothed_(tarEleIndex).faceCentres;
	iNodeCords = nodeCoordsSmoothed_(eNodMatSmoothed_(tarEleIndex,:),:);
	%%compute direction vectors from iCoord to face centers as reference vectors
	refVec = iElefaceCentres - iCoord; %% dir vecs from volume center to face centers
	refVec2Vertices = iNodeCords - iCoord; 
	refVecNorm = vecnorm(refVec,2,2);
	refVec2VerticesNorm = vecnorm(refVec2Vertices,2,2);
	if 6==numFacesPerEleSmoothed_
		thresholdAng = 91.0;
	else
		thresholdAng = 90.0;
	end
	if isempty(find(0==[refVecNorm; refVec2VerticesNorm])) %% iCoord does NOT coincides with a vertex or face center
		refVec = refVec ./ refVecNorm; 
		normVecs = eleStructSmoothed_(tarEleIndex).faceNormals;
		
		%%compute angle deviation
		angleDevs = zeros(numFacesPerEleSmoothed_,1);
		for ii=1:numFacesPerEleSmoothed_
			angleDevs(ii) = acos(refVec(ii,:)*normVecs(ii,:)')/pi*180;
		end
		%% iCoord is out of tarEleIndex, using the relaxed 91 instead of 90 for numerical instability
		maxAngle = max(angleDevs);
		if maxAngle > thresholdAng, opt = 0; end	
	end 
end

function val = ElementInterpolationInverseDistanceWeighting(coords, vtxEntity, ips)
	%% Inverse Distance Weighting
	%% coords --> element vertex coordinates, Matrix: [N-by-3] 
	%% vtxEntity --> entities on element vertics, Matrix: [N-by-M], e.g., M = 6 for 3D stress tensor
	%% ips --> to-be interpolated coordinate, Vector: [1-by-3]
	
	e = -2;
	D = vecnorm(ips-coords,2,2);
	[sortedD, sortedMapVec] = sort(D);
    if 0==sortedD(1)
        val = vtxEntity(sortedMapVec(1),:); return;
    end
	sortedVtxVals = vtxEntity(sortedMapVec,:);
	wV = sortedD.^e;
	V = sortedVtxVals.*wV;	
	val = sum(V) / sum(wV);
end

function principalStress = ComputePrincipalStress(cartesianStress)
	%% "cartesianStress" is in the order: Sigma_xx, Sigma_yy, Sigma_zz, Sigma_yz, Sigma_zx, Sigma_xy
	principalStress = zeros(1, 12);
	A = cartesianStress([1 6 5; 6 2 4; 5 4 3]);
	[eigenVec, eigenVal] = eig(A);
	principalStress([1 5 9]) = diag(eigenVal);
	principalStress([2 3 4 6 7 8 10 11 12]) = reshape(eigenVec,1,9);
end

function oPSLs = CompactStreamlines(iPSLs, truncatedThreshold)
	tarIndice = [];
	for ii=1:length(iPSLs)
		if iPSLs(ii).length > truncatedThreshold
			tarIndice(end+1,1) = ii;
		end
	end
	oPSLs = iPSLs(tarIndice);
	if isempty(oPSLs), oPSLs = []; end
end

function [nextElementIndex, p1, opt] = SearchNextIntegratingPointOnUnstructuredMesh(oldElementIndex, physicalCoordinates, sPoint, relocatingP1)
	global eleCentroidListSmoothed_; 
	global eNodMatSmoothed_; 
	global nodStructSmoothed_; 
	global eleStateSmoothed_;

	p1 = physicalCoordinates;
	nextElementIndex = oldElementIndex;	
	opt = IsThisPointWithinThatElement(oldElementIndex, p1);

	if opt
		return;
	else	
		tarNodes = eNodMatSmoothed_(oldElementIndex,:); 
		potentialElements = unique([nodStructSmoothed_(tarNodes(:)).adjacentEles]);
		adjEleCtrs = eleCentroidListSmoothed_(potentialElements,:);
		disList = vecnorm(p1-adjEleCtrs, 2, 2);
		[~, reSortMap] = sort(disList);
		potentialElements = potentialElements(reSortMap);
		for jj=1:length(potentialElements)
			iEle = potentialElements(jj);
			opt = IsThisPointWithinThatElement(iEle, p1);
			if opt, nextElementIndex = iEle; return; end
		end		
	end
	%%Scaling down the stepsize via Dichotomy
	if relocatingP1 && 0==opt	
		nn = 5;
		ii = 1;	
		while ii<=nn
			p1 = (sPoint+p1)/2;
			disList = vecnorm(p1-adjEleCtrs, 2, 2);
			[~, reSortMap] = sort(disList);
			potentialElements = potentialElements(reSortMap);			
			for jj=1:length(potentialElements)
				iEle = potentialElements(jj);
				opt = IsThisPointWithinThatElement(iEle, p1);
				if opt, nextElementIndex = iEle; return; end
			end
			ii = ii + 1;
		end
		if 0==eleStateSmoothed_(oldElementIndex)
			nextElementIndex = oldElementIndex; opt = 1;
		end
	end	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImportStressFields(fileName)
	global boundingBoxSmoothed_;
	global meshTypeSmoothed_;
	global numNodesSmoothed_;
	global nodeCoordsSmoothed_;
	global numElesSmoothed_;
	global eNodMatSmoothed_;
	global numNodesPerEleSmoothed_;
	global numEdgesPerEleSmoothed_;
	global eleStateSmoothed_;
	global nodStateSmoothed_;
	global boundaryEdgeNodMatSmoothed_;
	global cartesianStressFieldSmoothed_;
	global principalStressFieldSmoothed_;
	global loadingCondSmoothed_; 
	global fixingCondSmoothed_;	
	global eleCentroidListSmoothed_;
	global silhouetteStructSmoothed_;
	global eleSizeListSmoothed_;
	
	global nodStructSmoothed_; 
	global eleStructSmoothed_; 
	global boundaryElementsSmoothed_; 
	global boundaryNodesSmoothed_;
	global numFacesPerEleSmoothed_;
	global numNodesPerFaceSmoothed_;
	global eleFaceIndicesSmoothed_;	
	%%Read mesh and cartesian stress field
	fid = fopen(fileName, 'r');
	%%Mesh
	fgetl(fid); 
	domainType = fscanf(fid, '%s', 1);
	if ~strcmp(domainType, 'Solid'), warning('Un-supported Data!'); return; end
	meshTypeSmoothed_ = fscanf(fid, '%s', 1);
	if ~(strcmp(meshTypeSmoothed_, 'Hex') || strcmp(meshTypeSmoothed_, 'Tet')), warning('Un-supported Mesh!'); return; end
	meshOrder = fscanf(fid, '%d', 1);
	if 1~=meshOrder, warning('Un-supported Mesh!'); return; end
	startReadingVertices = fscanf(fid, '%s', 1);
	if ~strcmp(startReadingVertices, 'Vertices:'), warning('Un-supported Data!'); return; end
	numNodesSmoothed_ = fscanf(fid, '%d', 1);
	nodeCoordsSmoothed_ = fscanf(fid, '%e %e %e', [3, numNodesSmoothed_])'; 
	startReadingElements = fscanf(fid, '%s', 1);
	if ~strcmp(startReadingElements, 'Elements:'), warning('Un-supported Data!'); return; end
	numElesSmoothed_ = fscanf(fid, '%d', 1);
	switch meshTypeSmoothed_
		case 'Hex'
			eNodMatSmoothed_ = fscanf(fid, '%d %d %d %d %d %d %d %d', [8, numElesSmoothed_])'; 
		case 'Tet'
			eNodMatSmoothed_ = fscanf(fid, '%d %d %d %d', [4, numElesSmoothed_])'; 
	end

	startReadingLoads = fscanf(fid, '%s %s', 2); 
	if ~strcmp(startReadingLoads, 'NodeForces:'), warning('Un-supported Data!'); return; end
	numLoadedNodes = fscanf(fid, '%d', 1);
	if numLoadedNodes>0, loadingCondSmoothed_ = fscanf(fid, '%d %e %e %e', [4, numLoadedNodes])'; else, loadingCondSmoothed_ = []; end
    
	startReadingFixations = fscanf(fid, '%s %s', 2);
    if ~strcmp(startReadingFixations, 'FixedNodes:'), warning('Un-supported Data!'); return; end
	numFixedNodes = fscanf(fid, '%d', 1);
	if numFixedNodes>0, fixingCondSmoothed_ = fscanf(fid, '%d', [1, numFixedNodes])'; else, fixingCondSmoothed_ = []; end
    
	startReadingStress = fscanf(fid, '%s %s', 2); 
	if ~strcmp(startReadingStress, 'CartesianStress:'), warning('Un-supported Data!'); return; end
	numValidEles = fscanf(fid, '%d', 1);
	principalStressFieldSmoothed_ = fscanf(fid, '%e %e %e %e %e %e %e %e %e %e %e %e', [12, numValidEles])';		
	fclose(fid);
		
	%%Extract Boundary Element Info.
	switch meshTypeSmoothed_
		case 'Hex'
			numNodesPerEleSmoothed_ = 8;
			numFacesPerEleSmoothed_ = 6;
			numNodesPerFaceSmoothed_ = 4;
			eleFaceIndicesSmoothed_ = [4 3 2 1; 5 6 7 8; 1 2 6 5; 8 7 3 4; 5 8 4 1; 2 3 7 6];
			pp = [1 2 3 4];	
		case 'Tet'
			numNodesPerEleSmoothed_ = 4;
			numFacesPerEleSmoothed_ = 4;
			numNodesPerFaceSmoothed_ = 3;
			eleFaceIndicesSmoothed_ = [3 2 1;  1 2 4;  2 3 4;  3 1 4];
			pp = [1 2 2 3];	
	end
	
	boundingBoxSmoothed_ = [min(nodeCoordsSmoothed_, [], 1); max(nodeCoordsSmoothed_, [], 1)];	
	[surfaceMeshElements_, surfaceMeshNodeCoords_, nodStateSmoothed_, boundaryNodesSmoothed_] = ExtractBoundaryInfoFromSolidMesh();
	%%Extract Silhouette for Vis.
	silhouetteStructSmoothed_.vertices = surfaceMeshNodeCoords_;
	silhouetteStructSmoothed_.faces = surfaceMeshElements_;

	%%element centroids
	eleNodCoordListX = nodeCoordsSmoothed_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMatSmoothed_);
	eleNodCoordListY = nodeCoordsSmoothed_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMatSmoothed_);
	eleNodCoordListZ = nodeCoordsSmoothed_(:,3); eleNodCoordListZ = eleNodCoordListZ(eNodMatSmoothed_);
	eleCentroidListSmoothed_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2) sum(eleNodCoordListZ,2)]/numNodesPerEleSmoothed_;		
	
	%% Build Element Tree for Unstructured Quad-Mesh	
	iNodStruct = struct('adjacentEles', []); 
	nodStructSmoothed_ = repmat(iNodStruct, numNodesSmoothed_, 1);
	for ii=1:numElesSmoothed_
		for jj=1:numNodesPerEleSmoothed_
			nodStructSmoothed_(eNodMatSmoothed_(ii,jj)).adjacentEles(1,end+1) = ii;
		end
	end		
	boundaryElementsSmoothed_ = unique([nodStructSmoothed_(boundaryNodesSmoothed_).adjacentEles]);
	boundaryElementsSmoothed_ = boundaryElementsSmoothed_(:);
	eleStateSmoothed_ = zeros(numElesSmoothed_,1);
	eleStateSmoothed_(boundaryElementsSmoothed_,1) = 1;		
	%% build element tree		
	iEleStruct = struct('faceCentres', [], 'faceNormals', []); %%pure-Hex
	eleStructSmoothed_ = repmat(iEleStruct, numElesSmoothed_, 1);
	for ii=1:numElesSmoothed_
		iNodes = eNodMatSmoothed_(ii,:);
		iEleVertices = nodeCoordsSmoothed_(iNodes, :);
		iEleFacesX = iEleVertices(:,1); iEleFacesX = iEleFacesX(eleFaceIndicesSmoothed_);
		iEleFacesY = iEleVertices(:,2); iEleFacesY = iEleFacesY(eleFaceIndicesSmoothed_);
		iEleFacesZ = iEleVertices(:,3); iEleFacesZ = iEleFacesZ(eleFaceIndicesSmoothed_);
		ACs = [iEleFacesX(:,pp(1))-iEleFacesX(:,pp(3)) iEleFacesY(:,pp(1))-iEleFacesY(:,pp(3)) iEleFacesZ(:,pp(1))-iEleFacesZ(:,pp(3))];
		BDs = [iEleFacesX(:,pp(2))-iEleFacesX(:,pp(4)) iEleFacesY(:,pp(2))-iEleFacesY(:,pp(4)) iEleFacesZ(:,pp(2))-iEleFacesZ(:,pp(4))];
		iACxBD = cross(ACs,BDs); 
		aveNormal = iACxBD ./ vecnorm(iACxBD,2,2);			
		tmp = iEleStruct;			
		%% tmp.faceNormals = aveNormal;
		%% in case the node orderings on each element face are not constant
		tmp.faceCentres = [sum(iEleFacesX,2) sum(iEleFacesY,2) sum(iEleFacesZ,2)]/numNodesPerFaceSmoothed_;
		iEleCt = eleCentroidListSmoothed_(ii,:);
		refVecs = iEleCt - tmp.faceCentres; refVecs = refVecs ./ vecnorm(refVecs,2,2);
		dirEval = acos(sum(refVecs .* aveNormal, 2));
		dirDes = ones(numFacesPerEleSmoothed_,1); dirDes(dirEval<pi/2) = -1;
		faceNormals = dirDes .* aveNormal;
		tmp.faceNormals = faceNormals;
		eleStructSmoothed_(ii) = tmp;
	end
	
	%% Evaluate Element Sizes
	tmpSizeList = zeros(numFacesPerEleSmoothed_, numElesSmoothed_);
	for ii=1:numElesSmoothed_
		tmpSizeList(:,ii) = vecnorm(eleCentroidListSmoothed_(ii,:)-eleStructSmoothed_(ii).faceCentres,2,2);
	end
	eleSizeListSmoothed_ = 2*min(tmpSizeList,[],1)';		
end

function [boundaryFaceNodMat, boundaryFaceNodeCoords, nodState, boundaryNodes] = ExtractBoundaryInfoFromSolidMesh()
	global meshTypeSmoothed_;
	global numNodesSmoothed_;
	global numElesSmoothed_;
	global eNodMatSmoothed_;
	global nodeCoordsSmoothed_;
	global numNodesPerEleSmoothed_;
	global numFacesPerEleSmoothed_;
	global numNodesPerFaceSmoothed_;
	global eleFaceIndicesSmoothed_;

	eleFaces = eleFaceIndicesSmoothed_'; eleFaces = eleFaces(:)';
	patchIndices = eNodMatSmoothed_(:,eleFaces)';
	patchIndices = reshape(patchIndices(:), numNodesPerFaceSmoothed_, numFacesPerEleSmoothed_*numElesSmoothed_)';	
	tmp = sort(patchIndices,2);
	[uniqueFaces, ia, ic] = unique(tmp, 'stable', 'rows');
	leftFaceIDs = (1:numFacesPerEleSmoothed_*numElesSmoothed_)'; leftFaceIDs = setdiff(leftFaceIDs, ia);
	leftFaces = tmp(leftFaceIDs,:);
	[surfFaces, surfFacesIDsInUniqueFaces] = setdiff(uniqueFaces, leftFaces, 'rows');
	boundaryFaceNodMat = patchIndices(ia(surfFacesIDsInUniqueFaces),:);
	boundaryNodes = int32(unique(boundaryFaceNodMat));
	nodState = zeros(numNodesSmoothed_,1); nodState(boundaryNodes) = 1;	
	
	allNodes = zeros(numNodesSmoothed_,1);
	allNodes(boundaryNodes) = (1:numel(boundaryNodes))';
	boundaryFaceNodMat = allNodes(boundaryFaceNodMat);
	boundaryFaceNodeCoords = nodeCoordsSmoothed_(boundaryNodes,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowPSLs(varargin)
    global boundingBoxSmoothed_;
    global silhouetteStructSmoothed_;
	global majorPSLpoolSmoothed_;
	global mediumPSLpoolSmoothed_;
	global minorPSLpoolSmoothed_;
	
	numTarMajorPSLs = numel(majorPSLpoolSmoothed_);
	numTarMediumPSLs = numel(mediumPSLpoolSmoothed_);
	numTarMinorPSLs = numel(minorPSLpoolSmoothed_);
	color4MajorPSLs = struct('arr', []); color4MajorPSLs = repmat(color4MajorPSLs, numTarMajorPSLs, 1);
	color4MediumPSLs = struct('arr', []); color4MediumPSLs = repmat(color4MediumPSLs, numTarMediumPSLs, 1);
	color4MinorPSLs = struct('arr', []); color4MinorPSLs = repmat(color4MinorPSLs, numTarMinorPSLs, 1);	
	for ii=1:numTarMajorPSLs
		color4MajorPSLs(ii).arr = ones(1, majorPSLpoolSmoothed_(ii).length);
	end
	for ii=1:numTarMediumPSLs
		color4MediumPSLs(ii).arr = ones(1, mediumPSLpoolSmoothed_(ii).length);
	end			
	for ii=1:numTarMinorPSLs
		color4MinorPSLs(ii).arr = ones(1, minorPSLpoolSmoothed_(ii).length);
	end	
	if 0==nargin, thicknessScaling = 150; else, thicknessScaling = varargin{1}; end
	lineWidthTube = min(boundingBoxSmoothed_(2,:)-boundingBoxSmoothed_(1,:))/thicknessScaling;
	[gridXmajor, gridYmajor, gridZmajor, gridCmajor, ~] = ExpandPSLs2Tubes(majorPSLpoolSmoothed_, color4MajorPSLs, lineWidthTube);
	[gridXmedium, gridYmedium, gridZmedium, gridCmedium, ~] = ExpandPSLs2Tubes(mediumPSLpoolSmoothed_, color4MediumPSLs, lineWidthTube);
	[gridXminor, gridYminor, gridZminor, gridCminor, ~] = ExpandPSLs2Tubes(minorPSLpoolSmoothed_, color4MinorPSLs, lineWidthTube);
	handleMajorPSL = []; handleMediumPSL = []; handleMinorPSL = [];
	%%Show silhouette
	hSilo = patch(gca, silhouetteStructSmoothed_); hold(gca, 'on');	
	if ~isempty(gridXmajor)
		hold(gca, 'on'); 
		handleMajorPSL = surf(gca, gridXmajor, gridYmajor, gridZmajor, gridCmajor);
	end
	% if ~isempty(gridXmedium)
	% 	hold(gca, 'on');
	% 	handleMediumPSL = surf(gca, gridXmedium, gridYmedium, gridZmedium, gridCmedium);
	% end
	if ~isempty(gridXminor)
		hold(gca, 'on');
		handleMinorPSL = surf(gca, gridXminor, gridYminor, gridZminor, gridCminor);
	end	
	set(handleMajorPSL, 'FaceColor', [191 129 45]/255, 'EdgeColor', 'None'); %set(handleMajorPSL, 'FaceColor', [255 240 214]/255);
	%set(handleMediumPSL, 'FaceColor', [148 0 211]/255, 'EdgeColor', 'None'); %set(handleMediumPSL, 'FaceColor', [255 240 214]/255);
	set(handleMinorPSL, 'FaceColor', [53 151 143]/255, 'EdgeColor', 'None'); %set(handleMinorPSL, 'FaceColor', [255 240 214]/255);
	set(hSilo, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'None');
	view(gca, 30,0);
	camproj(gca, 'perspective');
	axis(gca, 'equal'); 
	axis(gca, 'tight');
	axis(gca, 'off');
	%%Lighting, Reflection
	lighting(gca, 'gouraud');
	material(gca, 'dull');
	camlight(gca, 'headlight', 'infinite');    
end

function [gridX, gridY, gridZ, gridC, gridIndices] = ExpandPSLs2Tubes(PSLs, colorSrc, r)
	%%Syntax
	%% [gridX, gridY, gridZ, gridC] = ExpandPSLs2Tubes(PSLs, colorSrc, r)
	gridX = [];
	gridY = [];
	gridZ = [];
	gridC = [];
	gridIndices = [];
	if isempty(PSLs), return; end
	n = 8; 
	numLines = length(PSLs);
	gridXYZ = zeros(3,n+1,1);
	gridC = zeros(n+1,1);
	gridIndices = struct('arr', []);
	gridIndices = repmat(gridIndices, numLines, 1);	
	for ii=1:numLines		
		curve = PSLs(ii).phyCoordList';
		npoints = size(curve,2);
		%deltavecs: average for internal points. first strecth for endpoitns.		
		dv = curve(:,[2:end,end])-curve(:,[1,1:end-1]);		
		%make nvec not parallel to dv(:,1)
		nvec=zeros(3,1); 
		[~,idx]=min(abs(dv(:,1))); 
		nvec(idx)=1;
		%precalculate cos and sing factors:
		cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
		sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
		%Main loop: propagate the normal (nvec) along the tube
		xyz = zeros(3,n+1,npoints+2);
		for k=1:npoints
			convec=cross(nvec,dv(:,k));
			convec=convec./norm(convec);
			nvec=cross(dv(:,k),convec);
			nvec=nvec./norm(nvec);
			%update xyz:
			xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1]) + cfact.*repmat(r*nvec,[1,n+1]) + sfact.*repmat(r*convec,[1,n+1]);
        end
		%finally, cap the ends:
		xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
		xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
		gridIndices(ii).arr = size(gridXYZ,3) : size(gridXYZ,3)+size(xyz,3)-1;
		gridXYZ(:,:,end+1:end+npoints+2) = xyz;	
		color = colorSrc(ii).arr;	
		c = [color(1) color color(end)];
		c = repmat(c, n+1, 1);
		gridC(:,end+1:end+npoints+2) = c;
	end		
	gridX = squeeze(gridXYZ(1,:,:)); 
	gridX(:,1) = [];
	gridY = squeeze(gridXYZ(2,:,:)); 
	gridY(:,1) = [];
	gridZ = squeeze(gridXYZ(3,:,:)); 
	gridZ(:,1) = [];
	gridC(:,1) = [];
end

function ShowProblemDescription()
	global nodeCoordsSmoothed_;
	global loadingCondSmoothed_;
	global fixingCondSmoothed_;
	global boundingBoxSmoothed_;
	global silhouetteStructSmoothed_;

	
	hd = patch(silhouetteStructSmoothed_); hold('on');
	set(hd, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 0.9, 'EdgeColor', 'k');
    if ~isempty(loadingCondSmoothed_)
		lB = 0.2; uB = 1.0;
		amps = vecnorm(loadingCondSmoothed_(:,2:end),2,2);
		maxAmp = max(amps); minAmp = min(amps);
		if abs(minAmp-maxAmp)/(minAmp+maxAmp)<0.1
			scalingFac = 1;
		else
			if minAmp/maxAmp>lB/uB, lB = minAmp/maxAmp; end
			scalingFac = lB + (uB-lB)*(amps-minAmp)/(maxAmp-minAmp);
		end
		loadingDirVec = loadingCondSmoothed_(:,2:end)./amps.*scalingFac;
		coordLoadedNodes = nodeCoordsSmoothed_(loadingCondSmoothed_(:,1),:);
		amplitudesF = mean(boundingBoxSmoothed_(2,:)-boundingBoxSmoothed_(1,:))/5 * loadingDirVec;
		hold('on'); quiver3(coordLoadedNodes(:,1), coordLoadedNodes(:,2), coordLoadedNodes(:,3), amplitudesF(:,1), ...
			amplitudesF(:,2), amplitudesF(:,3), 0, 'Color', [255 127 0.0]/255, 'LineWidth', 2, 'MaxHeadSize', 1, 'MaxHeadSize', 1);
	end
    if ~isempty(fixingCondSmoothed_)
		tarNodeCoord = nodeCoordsSmoothed_(fixingCondSmoothed_(:,1),:);
		hold('on'); hd1 = plot3(tarNodeCoord(:,1), tarNodeCoord(:,2), tarNodeCoord(:,3), 'x', ...
			'color', [153 153 153]/255, 'LineWidth', 3, 'MarkerSize', 15);		
	end
	view(gca, 3);
	camproj(gca, 'perspective');
	axis(gca, 'equal'); 
	axis(gca, 'tight');
	axis(gca, 'off');
	
	%%Lighting, Reflection
	lighting(gca, 'gouraud');
	material(gca, 'dull');
	camlight(gca, 'headlight', 'infinite');
end