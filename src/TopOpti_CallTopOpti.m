function TopOpti_CallTopOpti(axHandle)
	global meshHierarchy_;
	global F_;
	global U_;	
	global constraintType_;
	global startingGuess_;
	global V_;
	global alphaMin_;
	global passiveElements_;
	
	%%0.
	if isempty(F_), FEA_ApplyBoundaryCondition(); end
	if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
	U_ = zeros(meshHierarchy_(1).numDOFs,1);	
	
	%%1. Apply for Boundary Condition
	FEA_ApplyBoundaryCondition();
	
	%%2. Setup FEA
	FEA_SetupVoxelBased();
	
	%%3. Run Optimization
	tStart2 = tic;
	switch constraintType_
		case 'Global'
			startingGuess_ = repmat(V_, meshHierarchy_(1).numElements, 1);
			startingGuess_(passiveElements_) = 1;			
			TopOpti_GlobalVolumeConstraint(axHandle);
		case 'Local'
			startingGuess_ = repmat(alphaMin_, meshHierarchy_(1).numElements, 1);
			startingGuess_(passiveElements_) = 1;			
			TopOpti_LocalVolumeConstraint(axHandle);
	end
	disp(['Running TopOpti Costs: ' sprintf('%10.3g',toc(tStart2)) 's']);
	IO_ExportTopOptiResults();
	disp('..........Done with Topology Optimization!');
end