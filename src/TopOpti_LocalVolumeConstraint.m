function TopOpti_LocalVolumeConstraint(axHandle)
	%%1. initialize inputting arguments	
	global DEBUG_; global outPath_;
	global MEXfunc_;
	global meshHierarchy_;
	global passiveElements_;
	global startingGuess_;
	global modulus_;
	global complianceDesign_; 
	global volumeFractionDesign_; 
	global complianceSolid_;
	
	global p_;
	global nLoop_;
	global rMin_;
	global minChange_;
	global maxSharpness_;	
	global move_; global beta_;    				
	global pMax_;
	global penalty_;
	global penaltyIncrement_;
	global penaltyUpdateIterations_;
	
	global cHist_; cHist_ = [];
	global volHist_; volHist_ = [];
	global sharpHist_; sharpHist_ = [];
	global consHist_; consHist_ = [];
	global lssIts_; lssIts_ = [];
	global tHist_; tHist_ = [];
    global densityLayout_; densityLayout_ = [];
	
	numElements = meshHierarchy_(1).numElements;
	%%2. prepare filter, remove checkerboard patterns
	tDensityFilteringClock = tic;
	if MEXfunc_
		rMin = rMin_;
		eleMapForward = meshHierarchy_(1).eleMapForward;
		resX = meshHierarchy_(1).resX;
		resY = meshHierarchy_(1).resY;
		resZ = meshHierarchy_(1).resZ;
		Hs = TopOpti_SetupDensityFilter_mex(rMin, numElements, eleMapForward, resX, resY, resZ);
	else
		TopOpti_BuildDensityFilter_matrixFree();
	end
	timeDensityFiltering = toc(tDensityFilteringClock);
	disp(['Building Density Filter Costs: ' sprintf('%10.3g', timeDensityFiltering) 's']);

	%% prepare PDE filter, local volume computation
	tLocalVolumeConstraintClock = tic;
	PDEmatrixFree = 1;
	if PDEmatrixFree
		TopOpti_BuildPDEfilter_matrixFree();
	else
		TopOpti_BuildPDEfilter();
	end
	timeLocalVolumeConstraintFiltering = toc(tLocalVolumeConstraintClock);
	disp(['Building PDE Filter Costs: ' sprintf('%10.3g',timeLocalVolumeConstraintFiltering) 's']);
	
	%%3. prepare optimizer
	startingGuess_ = double(startingGuess_);
	passiveElements = passiveElements_;	
	activeEles = (1:int32(numElements))'; activeEles = setdiff(activeEles,passiveElements);
	volMaxList = Initialize4GradedPorosity();
	x = startingGuess_;
	xTilde = x;
	xPhys = TopOpti_DualizeDesignVariable(xTilde);
	densityLayout_ = xPhys(:);
	fileName = sprintf(strcat(outPath_, 'intermeidateDensityLayout-It-%d.nii'), 0);
	IO_ExportDesignInVolume_nii(fileName);	
	
	xold1 = x(activeEles);	
	xold2 = xold1;
	low = 0;
	upp = 0;
	loopbeta = 0; 
	loop = 0;
	change = 1.0;
	sharpness = 1.0;
	onesArrSingle = ones(numElements,1);

	%%4. Evaluate Compliance of Fully Solid Domain
	meshHierarchy_(1).eleModulus = repmat(modulus_, 1, numElements);
	tSolvingFEAssemblingClock = tic;
	Solving_AssembleFEAstencil();
	itSolvingFEAssembling = toc(tSolvingFEAssemblingClock);
	tSolvingFEAiterationClock = tic;
	lssIts_(end+1,1) = Solving_CG_GMGS('printP_OFF');
	itSolvingFEAiteration = toc(tSolvingFEAiterationClock);
	ceList = TopOpti_ComputeUnitCompliance();
	complianceSolid_ = meshHierarchy_(1).eleModulus*ceList;
	disp(['Compliance of Fully Solid Domain: ' sprintf('%16.6e',complianceSolid_)]);
	tHist_(end+1,:) = [itSolvingFEAssembling itSolvingFEAiteration 0 0 0 0];
	disp([' It.: ' sprintf('%4i',0) ' Assembling Time: ', sprintf('%4i',itSolvingFEAssembling) 's;', ' Solver Time: ', sprintf('%4i',itSolvingFEAiteration) 's.']);	

	%% Conduct Stress Analysis on Solid Domain and Extract the Dominant Stress Directions
	disp('Stress Analysis on Solid Domain ...');
	tStressAnalysis = tic;
    [cartesianStressField, ~] = FEA_StressAnalysis();
	vonMisesStressPerElement = FEA_ComputePerElementVonMisesStress(cartesianStressField);
	dominantDirSolid = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressField); clear cartesianStressField
	vonMisesVolume = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	vonMisesVolume(meshHierarchy_(1).eleMapBack,1) = vonMisesStressPerElement;
	vonMisesVolume = reshape(vonMisesVolume, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	densityLayout_ = ones(size(densityLayout_));
	IO_ExportDesignWithOneProperty_nii(vonMisesVolume, strcat(outPath_, 'ResultVolume_Solid_vonMises.nii')); 
	niftiwrite(dominantDirSolid, strcat(outPath_, 'dominantDirSolid.nii'));
	disp(['Done with Stress Analysis (inc. extracting dominant stress directions) after ', sprintf('%.f', toc(tStressAnalysis)), 's']);
	
	%%5. optimization
	while loop < nLoop_ && change > minChange_ && sharpness>maxSharpness_
		perIteCost = tic;
		loopbeta = loopbeta+1; loop = loop+1; 
		
		%%5.1 & 5.2 FEA, objective and sensitivity analysis
		meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(xPhys);
	    tSolvingFEAssemblingClock = tic;
		Solving_AssembleFEAstencil();
		itSolvingFEAssembling = toc(tSolvingFEAssemblingClock);
	    tSolvingFEAiterationClock = tic;
		lssIts_(end+1,1) = Solving_CG_GMGS('printP_OFF');
		itSolvingFEAiteration = toc(tSolvingFEAiterationClock);
		ceList = TopOpti_ComputeUnitCompliance();
		itimeSolvingFEA = itSolvingFEAssembling + itSolvingFEAiteration;	
		
		tOptimizationClock = tic;
		ceNorm = ceList / complianceSolid_;
		cObj = meshHierarchy_(1).eleModulus * ceNorm;
		complianceDesign_ = cObj*complianceSolid_;
		volumeFractionDesign_ = double(sum(xPhys(:)) / numElements);
		dc = -TopOpti_DeMaterialInterpolation(xPhys).*ceNorm;
		itimeOptimization = toc(tOptimizationClock);	
		
		tLocalVolumeConstraintClock = tic;
		if PDEmatrixFree
			x_pde_hat = TopOpti_PDEFiltering_matrixFree(xPhys);
		else
			x_pde_hat = TopOpti_PDEFiltering(xPhys);
		end
		dfdx_pde = (sum(x_pde_hat.^p_ ./ volMaxList.^p_)/numElements)^(1/p_-1)*(x_pde_hat.^(p_-1) ./ volMaxList.^p_)/numElements;	
		itimeLocalVolumeConstraint = toc(tLocalVolumeConstraintClock);
		
		%%5.3 filtering/modification of sensitivity
		tDensityFilteringClock = tic;
		dx = TopOpti_DeDualizeDesignVariable(xTilde);
		if MEXfunc_
			dc = TopOpti_PerformDensityFiltering_mex(dc.*dx./Hs, rMin, numElements, eleMapForward, resX, resY, resZ);
		else
			dc = TopOpti_DensityFiltering_matrixFree(dc.*dx, 1);
		end		
		itimeDensityFiltering = toc(tDensityFilteringClock);
		
		%%5.4 solve the optimization probelm
		tLocalVolumeConstraintClock = tic;
		if PDEmatrixFree
			dfdx = TopOpti_PDEFiltering_matrixFree(dfdx_pde(:));
		else
			dfdx = TopOpti_PDEFiltering(dfdx_pde(:));
		end
		itimeLocalVolumeConstraint = itimeLocalVolumeConstraint + toc(tLocalVolumeConstraintClock);
		tDensityFilteringClock = tic;
		if MEXfunc_
			dfdx = TopOpti_PerformDensityFiltering_mex(dfdx(:).*dx./Hs, rMin, numElements, eleMapForward, resX, resY, resZ); dfdx = dfdx(:)';
		else
			dfdx = TopOpti_DensityFiltering_matrixFree(dfdx(:).*dx, 1)';
		end		
		itimeDensityFiltering = itimeDensityFiltering + toc(tDensityFilteringClock);
		
		tOptimizationClock = tic;
		m = 1;
		n = numel(activeEles);
		df0dx = dc;
		fval = (sum(x_pde_hat.^p_ ./ volMaxList.^p_)/numElements)^(1/p_) - 1;	
		xval_MMA = x(activeEles);
		df0dx_MMA = df0dx(activeEles);
		dfdx_MMA = dfdx(:,activeEles);
		xmin_MMA = max(0.0,xval_MMA-move_);
		xmax_MMA = min(1,xval_MMA+move_);		
		if MEXfunc_
			[xmma_MMA, xold1, xold2] = MMA_mex(m, n, xval_MMA, xmin_MMA, xmax_MMA, xold1, xold2, df0dx_MMA, fval, dfdx_MMA(:));
		else
			[xmma_MMA, xold1, xold2] = MMAseq(m, n, xval_MMA, xmin_MMA, xmax_MMA, xold1, xold2, df0dx_MMA, fval, dfdx_MMA(:));
		end
		change = max(abs(xmma_MMA(:)-xval_MMA(:)));	
		x = onesArrSingle; x(activeEles) = xmma_MMA;		
		itimeOptimization = itimeOptimization + toc(tOptimizationClock);
		
		tDensityFilteringClock = tic;
		if MEXfunc_
			xTilde = TopOpti_PerformDensityFiltering_mex(x, rMin, numElements, eleMapForward, resX, resY, resZ); xTilde = xTilde./Hs;
		else
			xTilde = TopOpti_DensityFiltering_matrixFree(x, 0);
		end		
		xPhys = TopOpti_DualizeDesignVariable(xTilde);
		itimeDensityFiltering = itimeDensityFiltering + toc(tDensityFilteringClock);	
		xPhys(passiveElements) = 1;	
		sharpness = 4*sum(sum(xPhys.*(ones(numElements,1)-xPhys)))/numElements;
		densityLayout_ = xPhys(:);

    	if 1==loop || 0==mod(loop, 5)
		    fileName = sprintf(strcat(outPath_, 'intermeidateDensityLayout-It-%d.nii'), loop);
		    IO_ExportDesignInVolume_nii(fileName);
        end
		if DEBUG_		
			[az, el] = view(axHandle);
			cla(axHandle);
			Vis_ShowDesignByDensityLayoutInIsosurface(axHandle);
			view(axHandle, az, el);
			Vis_UserLighting(axHandle);
			drawnow;
		end
		itimeTotal = toc(perIteCost);

		%%5.5 write opti. history
		cHist_(loop,1) = complianceDesign_;
		volHist_(loop,1) = volumeFractionDesign_;
		consHist_(loop,:) = fval;
		sharpHist_(loop,1) = sharpness;
		iTimeStatistics = [itSolvingFEAssembling itSolvingFEAiteration itimeOptimization itimeDensityFiltering itimeLocalVolumeConstraint itimeTotal];
		tHist_(end+1,:) = iTimeStatistics;
		
		%%5.6 print results
		disp([' It.: ' sprintf('%i',loop) '... Obj.: ' sprintf('%10.4e',complianceDesign_) '; Vol.: ' sprintf('%6.3f',volumeFractionDesign_) ...
			 '; Sharp: ' sprintf('%10.4e',sharpness) '; Change: ' sprintf('%10.4e',change) '; Cons.: ' sprintf('%10.4e',fval)]);
		disp([' It.: ' sprintf('%i',loop) ' (Time)... Total per-It.: ' sprintf('%8.2e',itimeTotal), 's; Assemb.: ', ...
			sprintf('%8.2e',itSolvingFEAssembling), 's; CG: ', sprintf('%8.2e',itSolvingFEAiteration), ...
				's; Opti.: ', sprintf('%8.2e',itimeOptimization), 's; Filtering: ', sprintf('%8.2e',itimeDensityFiltering), ...
					's; LVF: ', sprintf('%8.2e',itimeLocalVolumeConstraint) 's.']);
			
		%%5.7 update Heaviside regularization parameter
		if beta_ < pMax_ && loopbeta >= 40
			beta_ = 2*beta_;
			loopbeta = 0;
			change = 1.0;
			sharpness = 1.0;
			fprintf('Parameter beta increased to %g.\n',beta_);			
		end

		%%5.8 update penalty for SIMP
		if ~mod(loop, penaltyUpdateIterations_)
			penalty_ = penalty_ + penaltyIncrement_;
			penalty_ = min([penalty_, 3.0]);
			fprintf('Penalty in SIMP increased to %g.\n', penalty_);	
		end			
	end

	disp('Stress Analysis on Design ...');
	tStressAnalysis = tic;
    [cartesianStressField, ~] = FEA_StressAnalysis();
	vonMisesStressPerElement = FEA_ComputePerElementVonMisesStress(cartesianStressField);
	dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressField); clear cartesianStressField
	niftiwrite(dominantDirDesign, strcat(outPath_, 'dominantDirDesign.nii'));
	vonMisesVolume = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	vonMisesVolume(meshHierarchy_(1).eleMapBack,1) = vonMisesStressPerElement;
	vonMisesVolume = reshape(vonMisesVolume, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	IO_ExportDesignWithOneProperty_nii(vonMisesVolume, strcat(outPath_, 'ResultVolume_Design_vonMises.nii')); 	
	disp(['Done with Stress Analysis (inc. extracting dominant stress directions) after ', sprintf('%.f', toc(tStressAnalysis)), 's']);

	disp('Compute Stress Aligment Scale between Solid and Design...');
	tStressAligmentAna = tic;
	alignmentMetricVolume = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
	IO_ExportDesignWithOneProperty_nii(alignmentMetricVolume, strcat(outPath_, 'ResultVolume_Design_StressAlignment.nii'));
	disp(['Done with Stress Alignment Analysis after ', sprintf('%.f', toc(tStressAligmentAna)), 's.']);
	
	fileName = strcat(outPath_, 'ResultVolume_Design.nii');
	IO_ExportDesignInVolume_nii(fileName);
	disp(['..........Solving FEA costs: ', sprintf('%10.4e', sum(sum(tHist_(:,1:2)))), 's.']);
	disp(['..........Optimization (inc. sentivity analysis, update) costs: ', sprintf('%10.4e', sum(tHist_(:,3))), 's.']);
	disp(['..........Performing Density-based Filtering costs: ', sprintf('%10.4e', sum(tHist_(:,4))), 's.']);
	disp(['..........Applying for Local Volume Constraint costs: ', sprintf('%10.4e', sum(tHist_(:,5))), 's.']);	
end

function gradedPorosityCtrlList = Initialize4GradedPorosity()
	global meshHierarchy_;
	global alphaMin_; global alphaMax_;
	global passiveElements_;
	
	alphaMax_ = alphaMin_;
	stepSize = (alphaMax_-alphaMin_) / (meshHierarchy_(1).resX-1);
	gradedPorosityCtrlList = alphaMin_*ones(meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	for ii=2:meshHierarchy_(1).resX					
		gradedPorosityCtrlList(:,ii,:) = gradedPorosityCtrlList(:,ii,:) + (ii-1)*stepSize;
	end	
	gradedPorosityCtrlList = gradedPorosityCtrlList(:);
	gradedPorosityCtrlList = gradedPorosityCtrlList(meshHierarchy_(1).eleMapBack);
	gradedPorosityCtrlList(passiveElements_) = 1;
end