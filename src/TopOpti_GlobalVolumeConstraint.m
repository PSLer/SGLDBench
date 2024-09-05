function TopOpti_GlobalVolumeConstraint(axHandle)
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
	
	global V_; 
	global nLoop_;
	global rMin_;
	global optimizer_;
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
	global tHist_; tHist_ = [];
	global lssIts_; lssIts_ = [];
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
	
	%%3. prepare optimizer
	passiveElements = passiveElements_;	
	activeEles = (1:int32(numElements))'; 
	activeEles = setdiff(activeEles,passiveElements);
	% x = single(startingGuess_);
	x = startingGuess_;
	xTilde = x;
	xPhys = TopOpti_DualizeDesignVariable(xTilde);
	densityLayout_ = xPhys(:);
	fileName = sprintf(strcat(outPath_, 'intermeidateDensityLayout-It-%d.nii'), 0);
	IO_ExportDesignInVolume_nii(fileName);	
	
	xold1 = x;	
	xold2 = x;
	low = 0;
	upp = 0;
	loopbeta = 0; 
	loop = 0;
	change = 1.0;
	sharpness = 1.0;
	% onesArrSingle = ones(numElements,1,'single');
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
	tHist_(end+1,:) = [itSolvingFEAssembling itSolvingFEAiteration 0 0 itSolvingFEAssembling + itSolvingFEAiteration];
	disp([' It.: ' sprintf('%4i',0) ' Assembling Time: ', sprintf('%4i',itSolvingFEAssembling) 's;', ' Solver Time: ', sprintf('%4i',itSolvingFEAiteration) 's.']);	
	
	%% Conduct Stress Analysis on Solid Domain and Extract the Dominant Stress Directions
	% disp('Stress Analysis on Solid Domain ...');
	% tStressAnalysis = tic;
	% dominantDirSolid = Common_ExtractDominantDirectionsFromPrincipalStressDirections();
	% niftiwrite(dominantDirSolid, strcat(outPath_, 'dominantDirSolid.nii'));
	% disp(['Done with Stress Analysis after ', sprintf('%.f', toc(tStressAnalysis)), 's']);
	
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
		
		tOptimizationClock = tic;
		ceNorm = ceList / complianceSolid_;
		cObj = meshHierarchy_(1).eleModulus * ceNorm;
		complianceDesign_ = cObj*complianceSolid_;
		volumeFractionDesign_ = double(sum(xPhys(:)) / numElements);
		dc = -TopOpti_DeMaterialInterpolation(xPhys).*ceNorm;
		dv = ones(numElements,1);
		itimeOptimization = toc(tOptimizationClock);
		
		%%5.3 filtering/modification of sensitivity
		tDensityFilteringClock = tic;
		dx = TopOpti_DeDualizeDesignVariable(xTilde);
		if MEXfunc_
			dc = TopOpti_PerformDensityFiltering_mex(dc.*dx./Hs, rMin, numElements, eleMapForward, resX, resY, resZ);
			dv = TopOpti_PerformDensityFiltering_mex(dv.*dx./Hs, rMin, numElements, eleMapForward, resX, resY, resZ);
		else
			dc = TopOpti_DensityFiltering_matrixFree(dc.*dx, 1);
			dv = TopOpti_DensityFiltering_matrixFree(dv.*dx, 1);			
		end
		itimeDensityFiltering = toc(tDensityFilteringClock);
		
		tOptimizationClock = tic;
		%%5.4 solve the optimization probelm
		switch optimizer_
			case 'MMA'
				m = 1;
				df0dx = dc;
				dfdx = dv'/(numElements*V_);
				iter = loopbeta;
				f0val = cObj;
				fval = sum(sum(xPhys)) / (numElements*V_) - 1;			
				
				a0 = 1;
				a = zeros(m,1);     
				c_ = ones(m,1)*1000;
				d = zeros(m,1);
				xval_MMA = x(activeEles);
				xold1_MMA = xold1(activeEles);
				xold2_MMA = xold2(activeEles);
				df0dx_MMA = df0dx(activeEles);
				df0dx2_MMA = zeros(numel(activeEles),1);
				dfdx_MMA = dfdx(:,activeEles');
				dfdx2_MMA = df0dx2_MMA';
				n = numel(activeEles);
				xmin_MMA = max(0.0,xval_MMA-move_);
				xmax_MMA = min(1,xval_MMA+move_);

				[xmma_MMA,~,~,~,~,~,~,~,~,low,upp] = ...
					mmasub(m,n,iter,xval_MMA,xmin_MMA,xmax_MMA,xold1_MMA,xold2_MMA,...
						f0val,df0dx_MMA,df0dx2_MMA,fval,dfdx_MMA,dfdx2_MMA,low,upp,a0,a,c_,d);				
				x = onesArrSingle;
				x(activeEles) = xmma_MMA;					
				xval = onesArrSingle; xval(activeEles) = xval_MMA;
				xold2 = xold1;
				xold1 = xval;
				change = max(x(:)-xval(:));
			case 'OC'		
				l1 = 0; l2 = 1e5;
				while (l2-l1)/(l1+l2) > 1e-3
					lmid = 0.5*(l2+l1);
					xnew = max(0,max(x-move_,min(1,min(x+move_,x.*sqrt(-dc./dv/lmid)))));
					xnew(passiveElements,1) = 1.0;
					fval = sum(xnew,1)/numElements;
					if fval > V_, l1 = lmid; else l2 = lmid; end				
				end		
				change = max(abs(xnew(:)-x(:)));
				x = xnew;	
		end
		itimeOptimization = itimeOptimization + toc(tOptimizationClock);

		tDensityFilteringClock = tic;
		if MEXfunc_
			xTilde = TopOpti_PerformDensityFiltering_mex(x, rMin, numElements, eleMapForward, resX, resY, resZ); xTilde = xTilde./Hs;
		else
			xTilde = TopOpti_DensityFiltering_matrixFree(x, 0);
		end
		itimeDensityFiltering = itimeDensityFiltering + toc(tDensityFilteringClock);
		xPhys = TopOpti_DualizeDesignVariable(xTilde);
		xPhys(passiveElements) = 1;		
		sharpness = 4*sum(sum(xPhys.*(ones(numElements,1)-xPhys)))/numElements;
		densityLayout_ = xPhys(:);

        %%if 1==loop || 0==mod(loop,5)
		    fileName = sprintf(strcat(outPath_, 'intermeidateDensityLayout-It-%d.nii'), loop);
		    IO_ExportDesignInVolume_nii(fileName);
        %%end
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
		iTimeStatistics = [itSolvingFEAssembling itSolvingFEAiteration itimeOptimization itimeDensityFiltering itimeTotal];
		tHist_(loop,:) = iTimeStatistics;
		
		%%5.6 print results
		disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4e',complianceDesign_) ' Vol.: ' sprintf('%6.3f',volumeFractionDesign_) ...
			 ' Sharp: ' sprintf('%10.4e',sharpness) ' Change: ' sprintf('%10.4e',change) ' Cons.: ' sprintf('%10.4e',fval)]);
		disp([' It.: ' sprintf('%4i',loop) ' Total Time per-It.: ' sprintf('%10.4e',itimeTotal) 's;', ' Assembling Time: ', ...
			sprintf('%10.4e',itSolvingFEAssembling) 's;', ' Solver Time: ', sprintf('%10.4e',itSolvingFEAiteration) 's;', ...
			' Optimization Time: ', sprintf('%10.4e',itimeOptimization) 's;', ' Filtering Time: ', sprintf('%10.4e',itimeDensityFiltering) 's.']);				 
		%%5.7 update Heaviside regularization parameter
		if beta_ < pMax_ && (loopbeta >= 40 || change <= 0.001)
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

	% disp('Stress Analysis on Design ...');
	% tStressAnalysis = tic;
	% dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections();
	% niftiwrite(dominantDirDesign, strcat(outPath_, 'dominantDirDesign.nii'));
	% disp(['Done with Stress Analysis after ', sprintf('%.f', toc(tStressAnalysis)), 's.']);
	
	% disp('Compute Stress Aligment Scale between Solid and Design...');
	% tStressAligmentAna = tic;
	% alignmentMetricVolume = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
	% niftiwrite(alignmentMetricVolume, strcat(outPath_, 'alignmentMetricVolume_byStress.nii'));
	% disp(['Done with Stress Alignment Analysis after ', sprintf('%.f', toc(tStressAligmentAna)), 's.']);
	
	fileName = strcat(outPath_, 'DesignVolume.nii');
	IO_ExportDesignInVolume_nii(fileName);
	disp(['..........Solving FEA costs: ', sprintf('%10.4e', sum(sum(tHist_(:,1:2)))), 's.']);
	disp(['..........Optimization (inc. sentivity analysis, update) costs: ', sprintf('%10.4e', sum(tHist_(:,3))), 's.']);
	disp(['..........Performing Density-based Filtering costs: ', sprintf('%10.4e', sum(tHist_(:,4))), 's.']);
end