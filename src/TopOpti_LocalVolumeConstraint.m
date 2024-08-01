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
	global move_; global beta_;    				
	global pMax_;
	
	global cHist_; cHist_ = [];
	global volHist_; volHist_ = [];
	global sharpHist_; sharpHist_ = [];
	global consHist_; consHist_ = [];
	global densityLayout_; densityLayout_ = [];
	global densityLayoutWithoutBoundary_; densityLayoutWithoutBoundary_ = [];
	
	%%2. prepare filter, remove checkerboard patterns
	tStart1 = tic;
	TopOpti_BuildDensityFilter();
	disp(['Building Density Filter Costs: ' sprintf('%10.3g',toc(tStart1)) 's']);

	%% prepare PDE filter, local volume computation
	tStart1 = tic;
	TopOpti_BuildPDEfilter();
	disp(['Building PDE Filter Costs: ' sprintf('%10.3g',toc(tStart1)) 's']);
	
	%%3. prepare optimizer
	numElements = meshHierarchy_(1).numElements;
	activeEles = (1:numElements)'; activeEles = setdiff(activeEles,passiveElements_);
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
	change = 1;

	%%4. Evaluate Compliance of Fully Solid Domain
	meshHierarchy_(1).eleModulus = repmat(modulus_, 1, numElements);
	ceList = TopOpti_ComputeUnitCompliance('printP_OFF');
	complianceSolid_ = meshHierarchy_(1).eleModulus*ceList;
	disp(['Compliance of Fully Solid Domain: ' sprintf('%16.6e',complianceSolid_)]);

	%%5. optimization
	while loop < nLoop_
		loopbeta = loopbeta+1; loop = loop+1; 
		
		%%5.1 & 5.2 FEA, objective and sensitivity analysis
		meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(xPhys);
		ceList = TopOpti_ComputeUnitCompliance('printP_OFF');
		ceNorm = ceList / complianceSolid_;
		cObj = meshHierarchy_(1).eleModulus * ceNorm;
		complianceDesign_ = cObj*complianceSolid_;
		volumeFractionDesign_ = sum(xPhys(:)) / numElements;
		dc = -TopOpti_DeMaterialInterpolation(xPhys).*ceNorm;
		x_pde_hat = TopOpti_PDEFiltering(xPhys);
		dfdx_pde = (sum(x_pde_hat.^p_ ./ volMaxList.^p_)/numElements)^(1/p_-1)*(x_pde_hat.^(p_-1) ./ volMaxList.^p_)/numElements;	

		%%5.3 filtering/modification of sensitivity
		dx = TopOpti_DeDualizeDesignVariable(xTilde);
		dc = TopOpti_DensityFiltering(dc.*dx, 1);
		
		%%5.4 solve the optimization probelm
		m = 1;
		df0dx = dc;
		df0dx2 = 0*df0dx;
		dfdx = TopOpti_PDEFiltering(dfdx_pde(:))';
		dfdx = TopOpti_DensityFiltering(dfdx'.*dx, 1)';
		dfdx2 = 0*dfdx;
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
		df0dx2_MMA = df0dx2(activeEles);
		dfdx_MMA = dfdx(:,activeEles');
		dfdx2_MMA = dfdx2(:,activeEles');
		n = length(activeEles);
		xmin_MMA = max(0.0,xval_MMA-move_);
		xmax_MMA = min(1,xval_MMA+move_);
		[xmma_MMA,ymma,zmma,lam,xsi,tmpEta,mu,zet,s,low,upp] = ...
			mmasub(m,n,iter,xval_MMA,xmin_MMA,xmax_MMA,xold1_MMA,xold2_MMA,...
				f0val,df0dx_MMA,df0dx2_MMA,fval,dfdx_MMA,dfdx2_MMA,low,upp,a0,a,c_,d);
		
		xmma = ones(numElements,1);
		xmma(activeEles) = xmma_MMA;
		xold1 = ones(numElements,1);
		xold1(activeEles) = xold1_MMA;
		xval = ones(numElements,1);
		xval(activeEles) = xval_MMA;
		xnew = xmma;
		xold2 = xold1;
		xold1 = xval;	
		
		xTilde = TopOpti_DensityFiltering(xnew, 0);
		xPhys = TopOpti_DualizeDesignVariable(xTilde);
		
		xPhys(passiveElements_) = 1;
		change = max(abs(xnew-x));
		x = xnew;	
		
		%%5.5 write opti. history
		cHist_(loop,1) = complianceDesign_;
		volHist_(loop,1) = volumeFractionDesign_;
		consHist_(loop,:) = fval;
		sharpHist_(loop,1) = 4*sum(sum(xPhys.*(ones(numElements,1)-xPhys)))/numElements;
		densityLayout_ = reshape(xPhys, numel(xPhys), 1);
		densityLayoutWithoutBoundary_ = densityLayout_;
		fileName = sprintf(strcat(outPath_, 'intermeidateDensityLayout-It-%d.mat'), loop);
		save(fileName, 'densityLayout_');	

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
			 ' Sharp: ' sprintf('%10.4e',sharpHist_(loop,1)) ' Cons.: ' sprintf('%10.4e',fval)]);

		%%5.7 update Heaviside regularization parameter
		if beta_ < pMax_ && loopbeta >= 40
			beta_ = 2*beta_;
			loopbeta = 0;
			change = 1;
			fprintf('Parameter beta increased to %g.\n',beta_);			
		end
	end
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