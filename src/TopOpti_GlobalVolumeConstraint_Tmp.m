function TopOpti_GlobalVolumeConstraint_Tmp()
	%% Aligning with [Traeff et al. CMAME, 2023]
	%%1. initialize inputting arguments	
	global outPath_;
	global meshHierarchy_;
	global startingGuess_;
	global complianceDesign_; 
	global volumeFractionDesign_; 
	global complianceSolid_;
	
	global V_; 
	global nLoop_;
	global rMin_;
	global minChange_;
	global maxSharpness_;
	global move_; 		
	
	global cHist_; cHist_ = [];
	global volHist_; volHist_ = [];
	global sharpHist_; sharpHist_ = [];
	global consHist_; consHist_ = [];
	global tHist_; tHist_ = [];
	global lssIts_; lssIts_ = [];
	global densityLayout_; densityLayout_ = [];
	
	numElements = meshHierarchy_(1).numElements;
	
	%%2. prepare filter, remove checkerboard patterns
	tDensityFilteringClock = tic;
	rMin = rMin_;
	eleMapForward = meshHierarchy_(1).eleMapForward;
	resX = meshHierarchy_(1).resX;
	resY = meshHierarchy_(1).resY;
	resZ = meshHierarchy_(1).resZ;
	Hs = TopOpti_SetupDensityFilter_mex(rMin, numElements, eleMapForward, resX, resY, resZ);
	timeDensityFiltering = toc(tDensityFilteringClock);
	disp(['Building Density Filter Costs: ' sprintf('%10.3g', timeDensityFiltering) 's']);
	
	%%3. prepare optimizer
	x = startingGuess_;
	xPhys = x;
	loop = 0;
	change = 1.0;
	sharpness = 1.0;
	complianceSolid_ = 1;
	%%5. optimization
	while loop < nLoop_ && change > minChange_ && sharpness>maxSharpness_
		perIteCost = tic;
		loop = loop+1; 
		
		%%5.1 & 5.2 FEA, objective and sensitivity analysis
		meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(xPhys);
		tSolvingFEAssemblingClock = tic;
	    Solving_AssembleFEAstencil();
		itSolvingFEAssembling = toc(tSolvingFEAssemblingClock);
		tSolvingFEAiterationClock = tic;
	    lssIts_(end+1,1) = Solving_CG_GMGS('printP_OFF');    
		itSolvingFEAiteration = toc(tSolvingFEAiterationClock);
		ceList = TopOpti_ComputeUnitCompliance();
		
		tOptimizationClock = tic;
		ceNorm = ceList / complianceSolid_;
		cObj = meshHierarchy_(1).eleModulus * ceNorm;
		complianceDesign_ = cObj*complianceSolid_;
		volumeFractionDesign_ = double(sum(xPhys(:)) / numElements);
		dc = -TopOpti_DeMaterialInterpolation(xPhys).*ceNorm;
		dv = ones(numElements,1);
		itimeOptimization = toc(tOptimizationClock);
		
		ft = 2;
		%%5.3 filtering/modification of sensitivity
		tDensityFilteringClock = tic;
		if 1==ft
			dc = TopOpti_PerformDensityFiltering_mex(x(:).*dc(:), rMin, numElements, eleMapForward, resX, resY, resZ)./Hs./max(1.0e-3,x(:));
		elseif 2==ft
			dc = TopOpti_PerformDensityFiltering_mex(dc./Hs, rMin, numElements, eleMapForward, resX, resY, resZ);
			dv = TopOpti_PerformDensityFiltering_mex(dv./Hs, rMin, numElements, eleMapForward, resX, resY, resZ);		
		end
		itimeDensityFiltering = toc(tDensityFilteringClock);
		
		tOptimizationClock = tic;
		%%5.4 solve the optimization probelm
		g = mean(xPhys(:))-V_;
		l1 = 0; l2 = 1e9;
		while (l2-l1)/(l1+l2) > 1e-6
			lmid = 0.5*(l2+l1);
			xnew = max(0,max(x-move_,min(1,min(x+move_,x.*sqrt(-dc./dv/lmid)))));
			gt=g+sum((dv(:).*(xnew(:)-x(:))));
			if gt>0, l1 = lmid; else l2 = lmid; end				
		end		
		change = max(abs(xnew(:)-x(:)));
		x = xnew;
		itimeOptimization = itimeOptimization + toc(tOptimizationClock);

		tDensityFilteringClock = tic;
		if 1==ft
			xPhys = xnew;
		elseif 2==ft
			xPhys = TopOpti_PerformDensityFiltering_mex(xnew, rMin, numElements, eleMapForward, resX, resY, resZ); xPhys = xPhys./Hs;
		end
		itimeDensityFiltering = itimeDensityFiltering + toc(tDensityFilteringClock);	
		sharpness = 4*sum(sum(xPhys.*(ones(numElements,1)-xPhys)))/numElements;
	
		itimeTotal = toc(perIteCost);
		
		%%5.5 write opti. history
		cHist_(loop,1) = complianceDesign_;
		volHist_(loop,1) = volumeFractionDesign_;
		consHist_(loop,:) = g;
		sharpHist_(loop,1) = sharpness;
		iTimeStatistics = [itSolvingFEAssembling itSolvingFEAiteration itimeOptimization itimeDensityFiltering itimeTotal];
		tHist_(loop,:) = iTimeStatistics;
		
		%%5.6 print results
		disp([' It.: ' sprintf('%i',loop) '... Obj.: ' sprintf('%10.4e',complianceDesign_) '; Vol.: ' sprintf('%6.3f',volumeFractionDesign_) ...
			 '; Sharp: ' sprintf('%10.4e',sharpness) ' Change: ' sprintf('%10.4e',change) '; Cons.: ' sprintf('%10.4e',g)]);
		disp([' It.: ' sprintf('%i',loop) ' (Time)... Total per-It.: ' sprintf('%8.2e',itimeTotal) 's;', ' Assemb.: ', ...
			sprintf('%8.2e',itSolvingFEAssembling), 's; CG: ', sprintf('%8.2e',itSolvingFEAiteration), ...
			's; Opti.: ', sprintf('%8.2e',itimeOptimization), 's; Filtering: ', sprintf('%8.2e',itimeDensityFiltering) 's.']);				 	
	end
	
	densityLayout_ = xPhys(:);
	fileName = strcat(outPath_, 'ResultVolume_Design.nii');
	IO_ExportDesignInVolume_nii(fileName);
	disp(['..........Solving FEA costs: ', sprintf('%10.4e', sum(sum(tHist_(:,1:2)))), 's.']);
	disp(['..........Optimization (inc. sentivity analysis, update) costs: ', sprintf('%10.4e', sum(tHist_(:,3))), 's.']);
	disp(['..........Performing Density-based Filtering costs: ', sprintf('%10.4e', sum(tHist_(:,4))), 's.']);
end