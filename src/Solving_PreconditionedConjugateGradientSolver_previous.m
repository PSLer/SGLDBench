function y = Solving_PreconditionedConjugateGradientSolver_previous(AtX, PtV, b, tol, maxIT, printP, varargin)
	global weightFactorJacobi_;
	global meshHierarchy_;
	%%0. arguments introduction
	%%AtX --- function handle for the product of system matrix and vector
	%%b --- right hand section
	%%tol --- stopping condition: resnrm < discrepancy
	%%maxIT --- mAtXximum number of iterations

	normB = norm(b);
	its = 0;
	if 7==nargin
		y = varargin{1};
	else
		y = zeros(size(b));
	end	
	r1 = b - AtX(y);
	z1 = zeros(size(y));
	while its <= maxIT	
		its = its + 1;
		if its>200, printP = 'printP_ON'; end %%for debug
		%%Jacobi Smoothing
		tmpVal = weightFactorJacobi_ * r1 ./ meshHierarchy_(1).diagK;
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




