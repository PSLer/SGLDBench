function tar = TopOpti_PDEFiltering_matrixFree(src)
	global TF_;
	global LF_;	
	global Permut_;		
	global KF_;  
	global maxIT_;
	global meshHierarchy_;
    global supMeshHierarchy4PDEfilter_;
	
	%%Element to Node
	tmpVal = zeros(meshHierarchy_(1).numNodes,1,'single');
if 1
	tStart = tic;
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,1)) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,1)) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,2)) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,2)) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,3)) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,3)) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,4)) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,4)) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,1)+1) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,1)+1) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,2)+1) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,2)+1) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,3)+1) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,3)+1) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,4)+1) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,4)+1) + src*(1/8);
	tEnd = toc(tStart)
else
	
end
    

	src = double(tmpVal);	
	
	%% Solving on Node
	%tar = Solving_PreconditionedConjugateGradientSolver4PDE(@Solving_KbyU_MatrixFree4PDE, @Solving_Vcycle4PDE, src, 1.0e-6, 200, 'printP_ON');
	tar = Solving_PreconditionedConjugateGradientSolver4PDE_nonPreCond(@Solving_KbyU_MatrixFree4PDE, src, 1.0e-6, 200, 'printP_ON');
	%%Node to Element
	tar = single(tar);
	tmpVal = zeros(meshHierarchy_(1).numElements,1,'single');
	blockIndex = Solving_MissionPartition(meshHierarchy_(1).numElements, 1.0e7);
	for jj=1:size(blockIndex,1)	
		rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))'; %%To avoid super-large data block
		iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
		iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
		tmpVal(rangeIndex,1) = sum(tar(iElesNodMat),2);
	end
	tar = tmpVal*(1/8);	
end

function y = Solving_PreconditionedConjugateGradientSolver4PDE_nonPreCond(AtX, b, tol, maxIT, printP)
	global weightFactorJacobi_;
	global supMeshHierarchy4PDEfilter_;
	%%0. arguments introduction
	%%AtX --- function handle for the product of system matrix and vector
	%%b --- right hand section
	%%tol --- stopping condition: resnrm < discrepancy
	%%maxIT --- mAtXximum number of iterations
	normB = norm(b);
	its = 0;
	y = zeros(size(b));	
	r1 = b - AtX(y);
	z1 = zeros(size(y));
	PtV = @(x) x;
	while its <= maxIT	
		its = its + 1;
		if its>200, printP = 'printP_ON'; end %%for debug
		
		z2 = PtV(r1);
		if 1==its
			p2 = z2;
		else
			beta = r1'*z2/(r0'*z1);
			p2 = z2 + beta*p1;			
		end
		tmpVal = AtX(p2);	
		
		alpha = r1'*z2/(p2'*tmpVal);
		y = y + alpha*p2;		
		r2 = r1 - alpha*tmpVal;
		
		resnorm = norm(r2)/normB;
		if strcmp(printP, 'printP_ON')
			disp([' It.: ' sprintf('%4i',its) ' Res.: ' sprintf('%16.6e',resnorm)]);
		end
		if resnorm<tol
			disp(['CG solver converged at iteration' sprintf('%5i', its) ' to a solution with relative residual' ...
					sprintf('%16.6e',resnorm)]);	
			break;
		end		
		%%update
		z1 = z2;
		p1 = p2;
		r0 = r1;
		r1 = r2;
	end	

	if its > maxIT
		warning('Exceed the maximum iterate numbers');
		disp(['The iterative process stops at residual = ' sprintf('%10.4f',resnorm)]);		
	end
end

function y = Solving_PreconditionedConjugateGradientSolver4PDE(AtX, PtV, b, tol, maxIT, printP)
	global weightFactorJacobi_;
	global supMeshHierarchy4PDEfilter_;
	%%0. arguments introduction
	%%AtX --- function handle for the product of system matrix and vector
	%%b --- right hand section
	%%tol --- stopping condition: resnrm < discrepancy
	%%maxIT --- mAtXximum number of iterations
	normB = norm(b);
	its = 0;
	y = zeros(size(b));	
	r1 = b - AtX(y);
	z1 = zeros(size(y));
	while its <= maxIT	
		its = its + 1;
		if its>200, printP = 'printP_ON'; end %%for debug
		%%Jacobi Smoothing
		tmpVal = weightFactorJacobi_ * r1 ./ supMeshHierarchy4PDEfilter_(1).diagK;
		y = y + tmpVal;
		r1 = r1 - AtX(tmpVal); clear tmpVal
		
		z2 = PtV(r1);
		if 1==its
			p2 = z2;
		else
			beta = r1'*z2/(r0'*z1);
			p2 = z2 + beta*p1;			
		end
		tmpVal = AtX(p2);	
		
		alpha = r1'*z2/(p2'*tmpVal);
		y = y + alpha*p2;		
		r2 = r1 - alpha*tmpVal;
		
		resnorm = norm(r2)/normB;
		if strcmp(printP, 'printP_ON')
			disp([' It.: ' sprintf('%4i',its) ' Res.: ' sprintf('%16.6e',resnorm)]);
		end
		if resnorm<tol
			disp(['CG solver converged at iteration' sprintf('%5i', its) ' to a solution with relative residual' ...
					sprintf('%16.6e',resnorm)]);	
			break;
		end		
		%%update
		z1 = z2;
		p1 = p2;
		r0 = r1;
		r1 = r2;
	end	

	if its > maxIT
		warning('Exceed the maximum iterate numbers');
		disp(['The iterative process stops at residual = ' sprintf('%10.4f',resnorm)]);		
	end
end

function productMV = Solving_KbyU_MatrixFree4PDE(uVec, varargin)	
	global meshHierarchy_;
	global supMeshHierarchy4PDEfilter_;
	if 1==nargin, iLevel = 1; else, iLevel = varargin{1}; end
	productMV = zeros(meshHierarchy_(iLevel).numNodes,1);
	Ks = supMeshHierarchy4PDEfilter_(iLevel).Ks;

	if 1==size(Ks,3)
		if 1~=iLevel, error('Wrong FEA Computing Stencil!'); end
		blockIndex = Solving_MissionPartition(meshHierarchy_(iLevel).numElements, 1.0e7);		
		for jj=1:size(blockIndex,1)	
			rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))'; %%To avoid super-large data block
			iElesNodMat = meshHierarchy_(iLevel).eNodMatHalf(rangeIndex,:);
			iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
			iIntermediateModulus = meshHierarchy_(iLevel).eleModulus(1,rangeIndex);
			subDisVec = uVec(iElesNodMat);
			subDisVec = subDisVec*Ks;
			productMV = productMV + accumarray(iElesNodMat(:),subDisVec(:),[meshHierarchy_(iLevel).numNodes 1]);
		end
	else
		eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(iLevel).eNodMatHalf);
		subDisVec = uVec(eNodMat);	
		for ii=1:meshHierarchy_(iLevel).numElements
			subDisVec(ii,:) = subDisVec(ii,:)*supMeshHierarchy4PDEfilter_(iLevel).Ks(:,:,ii);
		end
		productMV = productMV + accumarray(eNodMat(:),subDisVec(:),[meshHierarchy_(iLevel).numNodes 1]);
	end	
end

function deltaU = Solving_Vcycle4PDE(r)	
	global supMeshHierarchy4PDEfilter_;
	global weightFactorJacobi_;
	global numLevels_;
	global cholFacPDE_;
	global cholPermutPDE_;	
	
	%%0. preparation
	varVcycle = struct('x', 0, 'r', 0);
	varVcycle = repmat(varVcycle, numLevels_, 1);
	varVcycle(1).r = r;	

	%%1. Restriction. fine -> coarse
	for ii=2:numLevels_
		%%1.1 apply for smoother
		deltaX = weightFactorJacobi_ * varVcycle(ii-1).r ./ supMeshHierarchy4PDEfilter_(ii-1).diagK;
		varVcycle(ii-1).x = varVcycle(ii-1).x + deltaX;
		varVcycle(ii-1).r = varVcycle(ii-1).r - Solving_KbyU_MatrixFree4PDE(deltaX, ii-1); clear deltaX
		%%1.2 restrict residual
		varVcycle(ii).r = Solving_RestrictResidual4PDE(varVcycle(ii-1).r,ii);
	end	
	
	%%2. Directly solving on coarsest level
	varVcycle(end).x = (cholPermutPDE_*(cholFacPDE_'\(cholFacPDE_\(cholPermutPDE_'*varVcycle(end).r))));
	
	%%3. Interpolation. coarse -> fine
	for ii=numLevels_:-1:2
		%%3.1. interpolation			
		varVcycle(ii-1).x = varVcycle(ii-1).x + Solving_InterpolationDeviation4PDE(varVcycle(ii).x,ii);
		
		%%3.2 apply for smoother
		varVcycle(ii-1).x = varVcycle(ii-1).x + weightFactorJacobi_ * varVcycle(ii-1).r ./ supMeshHierarchy4PDEfilter_(ii-1).diagK;
	end
	deltaU = varVcycle(1).x;
end

function xFiner = Solving_InterpolationDeviation4PDE(xCoarser, ii)
	global meshHierarchy_;
	
	tmp = xCoarser(Common_RecoverHalfeNodMat(meshHierarchy_(ii).eNodMatHalf));
	tmp = meshHierarchy_(ii).multiGridOperatorRI * tmp';
	xFiner = accumarray(meshHierarchy_(ii).transferMat(:),tmp(:),[meshHierarchy_(ii).intermediateNumNodes 1]);
	xFiner = xFiner ./ meshHierarchy_(ii).transferMatCoeffi;
	xFiner = xFiner(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:);
end

function rCoaser = Solving_RestrictResidual4PDE(rFiner,ii)
	global meshHierarchy_;
	rFiner1 = zeros(meshHierarchy_(ii).intermediateNumNodes,1);
	rFiner1(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:) = rFiner;
	rFiner1 = rFiner1./meshHierarchy_(ii).transferMatCoeffi;
	eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(ii).eNodMatHalf);
	eNodMat = eNodMat(:);
	tmp = rFiner1(meshHierarchy_(ii).transferMat);
	tmp = tmp' * meshHierarchy_(ii).multiGridOperatorRI;
	rCoaser = accumarray(eNodMat,tmp(:),[meshHierarchy_(ii).numNodes 1]);
end