function FEA_ApplyBoundaryCondition()
	global meshHierarchy_;
	global F_;
	global U_;
	global loadingCond_;		
	global fixingCond_;
	
	%%Pre-Check
	[~, nodesLoadedFixed] = setdiff(fixingCond_(:,1), loadingCond_(:,1));
	fixingCond_ = fixingCond_(nodesLoadedFixed,:);
	[~,uniqueFixedNodes] = unique(fixingCond_(:,1));
	fixingCond_ = fixingCond_(uniqueFixedNodes,:);	
	[~,uniqueLoadedNodes] = unique(loadingCond_(:,1));
	loadingCond_ = loadingCond_(uniqueLoadedNodes,:);
	if isempty(loadingCond_), warning('No Loads!'); return; end
	if isempty(fixingCond_), warning('No Fixations!'); return; end
	
	%% Loading
	F_ = sparse(meshHierarchy_(1).numNodes, 3);
	F_(meshHierarchy_(1).nodesOnBoundary(loadingCond_(:,1)),:) = loadingCond_(:,2:end);
	F_ = reshape(F_',meshHierarchy_(1).numDOFs,1);
	U_ = zeros(meshHierarchy_(1).numDOFs,1);
	
	%%Fixing
	meshHierarchy_(1).fixedDOFs = [];
	meshHierarchy_(1).freeDOFs = [];
	fixedDOFs = 3*double(meshHierarchy_(1).nodesOnBoundary(fixingCond_(:,1)));
	fixedDOFs = fixedDOFs - [2 1 0];
	fixedDOFs = reshape(fixedDOFs', numel(fixedDOFs), 1);
	fixingState = fixingCond_(:,2:end)';
	fixedDOFs = fixedDOFs(1==fixingState(:)); %% E.g., X-dir is fixed, Y-dir is not
if 0	
	freeDOFs = setdiff((1:int32(meshHierarchy_(1).numDOFs))',fixedDOFs);
	meshHierarchy_(1).fixedDOFs = fixedDOFs;
	meshHierarchy_(1).freeDOFs = freeDOFs;
else
	freeDOFs = true(meshHierarchy_(1).numDOFs,1);
	freeDOFs(fixedDOFs) = false;
	meshHierarchy_(1).freeDOFs = freeDOFs;
	meshHierarchy_(1).fixedDOFs = false(meshHierarchy_(1).numDOFs,1);
	meshHierarchy_(1).fixedDOFs(fixedDOFs) = true;
end	
end
