function [y, varargout] = Solving_PreconditionedConjugateGradientSolver(AtX, PtV, b, tol, maxIT, printP, varargin)
	%%0. arguments introduction
	%%AtX --- function handle for the product of system matrix and vector
	%%b --- right hand section
	%%tol --- stopping condition: resnrm < discrepancy
	%%maxIT --- mAtXximum number of iterations
	normB = norm(b);
	if 7==nargin, y = varargin{1}; else, y = zeros(size(b)); end	
	rVec1 = b - AtX(y);	
	zVec = PtV(rVec1);	
	pVec = zVec;
	x1Val = zVec' * rVec1;
	for its=1:maxIT
		zVec = AtX(pVec);
		lambda = x1Val / (pVec' * zVec);
		y = y + lambda * pVec;
		rVec1 = rVec1 - lambda*zVec;
		resnorm = norm(rVec1)/normB;
		if strcmp(printP, 'printP_ON')
			disp([' It.: ' sprintf('%4i',its) ' Res.: ' sprintf('%16.6e',resnorm)]);
		end		
		if resnorm<tol
			disp(['CG solver converged at iteration' sprintf('%5i', its) ' to a solution with relative residual' sprintf('%16.6e',resnorm)]);
			break;
		end		
		zVec = PtV(rVec1);	
		x2Val = zVec' * rVec1;
		pVec = zVec + x2Val / x1Val * pVec;
		x1Val = x2Val;
	end

	if its == maxIT
		warning('Exceed the maximum iterate numbers');
		disp(['The iterative process stops at residual = ' sprintf('%10.4f',resnorm)]);		
	end
	if 2==nargout, varargout{1} = its; end
end
