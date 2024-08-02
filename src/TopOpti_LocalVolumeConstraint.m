function TopOpti_LocalVolumeConstraint(axHandle)
	%%1. initialize inputting arguments	
	global DEBUG_; global outPath_;
	global meshHierarchy_;
	global passiveElements_;
	global startingGuess_;
	global modulus_;
	global complianceDesign_; 
	global volumeFractionDesign_; 
	global complianceSolid_;
	
	global p_;
	global nLoop_;
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
	global densityLayout_; densityLayout_ = [];
	global densityLayoutWithoutBoundary_; densityLayoutWithoutBoundary_ = [];
	
	%%2. prepare filter, remove checkerboard patterns
	tStart1 = tic;
	TopOpti_BuildDensityFilter_matrixFree();
	disp(['Building Density Filter Costs: ' sprintf('%10.3g',toc(tStart1)) 's']);

	%% prepare PDE filter, local volume computation
	tStart1 = tic;
	TopOpti_BuildPDEfilter();
	disp(['Building PDE Filter Costs: ' sprintf('%10.3g',toc(tStart1)) 's']);
	
	%%3. prepare optimizer
	passiveElements = passiveElements_;
	numElements = meshHierarchy_(1).numElements;
	activeEles = (1:int32(numElements))'; activeEles = setdiff(activeEles,passiveElements);
	volMaxList = Initialize4GradedPorosity();
	x = startingGuess_;
	xTilde = x;
	xPhys = TopOpti_DualizeDesignVariable(xTilde);
	xold1 = x;	
	xold2 = x;
	low = 0;
	upp = 0;
	loopbeta = 0; 
	loop = 0;
	change = 1.0;
	sharpness = 1.0;
	onesArrSingle = ones(numElements,1,'single');

	%%4. Evaluate Compliance of Fully Solid Domain
	meshHierarchy_(1).eleModulus = repmat(modulus_, 1, numElements);
	ceList = TopOpti_ComputeUnitCompliance('printP_ON');
	complianceSolid_ = meshHierarchy_(1).eleModulus*ceList;
	disp(['Compliance of Fully Solid Domain: ' sprintf('%16.6e',complianceSolid_)]);

	%%5. optimization
	while loop < nLoop_ && change > minChange_ && sharpness>maxSharpness_
		loopbeta = loopbeta+1; loop = loop+1; 
		
		%%5.1 & 5.2 FEA, objective and sensitivity analysis
		meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(xPhys);
		ceList = TopOpti_ComputeUnitCompliance('printP_ON');
		ceNorm = ceList / complianceSolid_;
		cObj = meshHierarchy_(1).eleModulus * ceNorm;
		complianceDesign_ = cObj*complianceSolid_;
		volumeFractionDesign_ = double(sum(xPhys(:)) / numElements);
		dc = -TopOpti_DeMaterialInterpolation(xPhys).*ceNorm;
		x_pde_hat = TopOpti_PDEFiltering(xPhys);
		dfdx_pde = (sum(x_pde_hat.^p_ ./ volMaxList.^p_)/numElements)^(1/p_-1)*(x_pde_hat.^(p_-1) ./ volMaxList.^p_)/numElements;	

		%%5.3 filtering/modification of sensitivity
		dx = TopOpti_DeDualizeDesignVariable(xTilde);
		dc = TopOpti_DensityFiltering_matrixFree(dc.*dx, 1);
		
		%%5.4 solve the optimization probelm
		m = 1;
		df0dx = dc;
		dfdx = TopOpti_PDEFiltering(dfdx_pde(:))';
		dfdx = TopOpti_DensityFiltering_matrixFree(dfdx'.*dx, 1)';
		iter = loopbeta;
		f0val = cObj;	
		fval = (sum(x_pde_hat.^p_ ./ volMaxList.^p_)/numElements)^(1/p_) - 1;	
		
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
			mmasub(m,n,iter,double(xval_MMA),double(xmin_MMA),double(xmax_MMA),double(xold1_MMA),double(xold2_MMA),...
				f0val,double(df0dx_MMA),df0dx2_MMA,double(fval),double(dfdx_MMA),dfdx2_MMA,low,upp,a0,a,c_,d);
		
		x = onesArrSingle;
		x(activeEles) = single(xmma_MMA);	
		xval = onesArrSingle; xval(activeEles) = xval_MMA;
		xold2 = xold1;
		xold1 = xval;	
		change = max(abs(x(:)-xval(:)));
		
		xTilde = TopOpti_DensityFiltering_matrixFree(x, 0);
		xPhys = TopOpti_DualizeDesignVariable(xTilde);
		xPhys(passiveElements) = 1;	
		sharpness = 4*sum(sum(xPhys.*(ones(numElements,1)-xPhys)))/numElements;
		
		%%5.5 write opti. history
		cHist_(loop,1) = complianceDesign_;
		volHist_(loop,1) = volumeFractionDesign_;
		consHist_(loop,:) = fval;
		sharpHist_(loop,1) = sharpness;
		densityLayout_ = reshape(xPhys, numel(xPhys), 1);
		% fileName = sprintf(strcat(outPath_, 'intermeidateDensityLayout-It-%d.mat'), loop);
		% save(fileName, 'densityLayout_');
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
		
		%%5.6 print results
		disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4e',complianceDesign_) ' Vol.: ' sprintf('%6.3f',volumeFractionDesign_) ...
			 ' Sharp: ' sprintf('%10.4e',sharpness) ' Change: ' sprintf('%10.4e',change) ' Cons.: ' sprintf('%10.4e',fval)]);

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
	
	fileName = strcat(outPath_, 'DesignVolume.nii');
	IO_ExportDesignInVolume_nii(fileName);		
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