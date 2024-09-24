function z = Solving_Vcycle_NEW(r)	
	global meshHierarchy_;
	global weightFactorJacobi_;
	global numLevels_;
	global cholFac_; global cholPermut_;
	
	z = 0*r;
	z = smthdmpjac(z,r,l)
	
	
	
	
	
	
	
	
	
	%%0. preparation
	varVcycle = struct('x', 0, 'r', 0);
	varVcycle = repmat(varVcycle, numLevels_, 1);

	%%1. Restriction. fine -> coarse
	for ii=2:numLevels_
		if 2==ii
			%%1.1 apply for smoother
			rTilde = weightFactorJacobi_ * r ./ meshHierarchy_(ii-1).diagK;

			%%1.2 restrict residual
			varVcycle(ii).r = Solving_RestrictResidual(r,ii);		
		else
			%%1.1 apply for smoother
			varVcycle(ii-1).x = weightFactorJacobi_ * varVcycle(ii-1).r ./ meshHierarchy_(ii-1).diagK;
			
			%%1.2 restrict residual
			varVcycle(ii).r = Solving_RestrictResidual(varVcycle(ii-1).r,ii);		
		end
	end	
	
	%%2. Directly solving on coarsest level
	xCoarsest = zeros(meshHierarchy_(end).numDOFs,1);
	xCoarsest(meshHierarchy_(end).freeDOFs) = ...
		(cholPermut_*(cholFac_'\(cholFac_\(cholPermut_'*varVcycle(end).r(meshHierarchy_(end).freeDOFs)))));
	varVcycle(end).x = xCoarsest;
	
	%%3. Interpolation. coarse -> fine
	for ii=numLevels_:-1:2
		if 2==ii
			%%3.1. interpolation
			rTilde = rTilde + Solving_InterpolationDeviation(varVcycle(ii).x,ii);
			
			%%3.2 apply for smoother
			rTilde = rTilde + weightFactorJacobi_ * r ./ meshHierarchy_(ii-1).diagK;
		else
			%%3.1. interpolation			
			varVcycle(ii-1).x = varVcycle(ii-1).x + Solving_InterpolationDeviation(varVcycle(ii).x,ii);
			
			%%3.2 apply for smoother
			varVcycle(ii-1).x = varVcycle(ii-1).x + weightFactorJacobi_ * varVcycle(ii-1).r ./ meshHierarchy_(ii-1).diagK;		
		end
	end
	rTilde(meshHierarchy_(1).fixedDOFs) = 0;
	clear varVcycle
end

function u = Solving_JacobiSmoother_New(u, b, iLayer)
	global meshHierarchy_;
	global weightFactorJacobi_;
	
	omega = weightFactorJacobi_;
	invD = meshHierarchy_(iLayer).diagK;
	
	u = u - omega*(Solving_KbyU_MatrixFree(u, iLayer))./invD + omega*b./invD;
	
end

function [u] = smthdmpjac(u,A,b,invD,omega,nswp)
	% u = u - omega*invD.*(A*u) + omega*invD.*b;
	global meshHierarchy_;
	global weightFactorJacobi_;
	
	
	if 1==nargin, iLayer = 1; else, iLayer = varargin{1}; end
	deltaX = weightFactorJacobi_ * var.r ./ meshHierarchy_(iLayer).diagK;
	varargout{1} = var.x + deltaX;
	if 2==nargout
		varargout{2} = var.r - Solving_KbyU_MatrixFree_previous(deltaX, iLayer);
	end	
end