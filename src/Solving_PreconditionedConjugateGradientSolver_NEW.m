function [y, varargout] = Solving_PreconditionedConjugateGradientSolver_NEW(AtX, PtV, b, tol, maxIT, printP, varargin)
	%%0. arguments introduction
	%%AtX --- function handle for the product of system matrix and vector
	%%b --- right hand section
	%%tol --- stopping condition: resnrm < discrepancy
	%%maxIT --- mAtXximum number of iterations

	
	its = 0;
	if 7==nargin
		y = varargin{1};
	else
		y = zeros(size(b));
	end
	r = b - AtX(y);
	res0 = norm(b);

	% rTildeVec = PtV(r);

	% p = rTildeVec;

	while its <= maxIT	
		
		its = its + 1;
		z = PtV(r);
		rho = r' * z;
		
		if 1==its
			p = z;
		else
			beta=rho/rho_p;
			p=beta*p+z;
		end
		q = AtX(p);
		dpr = p'*q;
		alpha=rho/dpr;
		y=y+alpha*p;
		r=r-alpha*q;
		rho_p=rho;
		
		resnorm = norm(r)/res0;
		if strcmp(printP, 'printP_ON')
			disp([' It.: ' sprintf('%4i',its) ' Res.: ' sprintf('%16.6e',resnorm)]);
		end
		if resnorm<tol
			disp(['CG solver converged at iteration' sprintf('%5i', its) ' to a solution with relative residual' ...
					sprintf('%16.6e',resnorm)]);	
			break;
		end
	end	

	if its > maxIT
		warning('Exceed the maximum iterate numbers');
		disp(['The iterative process stops at residual = ' sprintf('%10.4f',resnorm)]);		
	end
	if 2==nargout, varargout{1} = its; end
	clear r rho_p p q
end

