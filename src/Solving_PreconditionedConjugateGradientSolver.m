function [y, varargout] = Solving_PreconditionedConjugateGradientSolver(AtX, PtV, b, tol, maxIT, printP, varargin)
	%%0. arguments introduction
	%%AtX --- function handle for the product of system matrix and vector
	%%b --- right hand section
	%%tol --- stopping condition: resnrm < discrepancy
	%%maxIT --- mAtXximum number of iterations
% global tMtV_; tMtV_ = 0;
% global tPtV_; tPtV_ = 0;
	normB = norm(b);
	its = 0;
	if 7==nargin
		y = varargin{1};
	else
		y = zeros(size(b));
	end
% tStart1 = tic;	
	rVec = b - AtX(y);
% tMtV_ = tMtV_ + toc(tStart1);
% tStart2 = tic;	
	rTildeVec = PtV(rVec);
% tPtV_ = tPtV_ + toc(tStart2);	
	pVec = rTildeVec;

	while its <= maxIT	
		its = its + 1;
% tStart1 = tic;		
		tmpVal = AtX(pVec);
% tMtV_ = tMtV_ + toc(tStart1);		
		% lambda = rTildeVec' * rVec / (pVec' * tmpVal);
		rTildeTimesrVec = rTildeVec' * rVec;
		lambda = rTildeTimesrVec / (pVec' * tmpVal);		
		y = y + lambda * pVec;
		r2Vec = rVec - lambda * tmpVal;
		resnorm = norm(r2Vec)/normB;
		if strcmp(printP, 'printP_ON')
			disp([' It.: ' sprintf('%4i',its) ' Res.: ' sprintf('%16.6e',resnorm)]);
		end
		if resnorm<tol
			disp(['CG solver converged at iteration' sprintf('%5i', its) ' to a solution with relative residual' ...
					sprintf('%16.6e',resnorm)]);	
			break;
		end
% tStart2 = tic;			
		r2TildeVec = PtV(r2Vec);
% tPtV_ = tPtV_ + toc(tStart2);
		% p2Vec = r2TildeVec + r2TildeVec' * r2Vec / (rTildeVec' * rVec) * pVec;
		p2Vec = r2TildeVec + r2TildeVec' * r2Vec / rTildeTimesrVec * pVec;
		%%update
		pVec = p2Vec;
		rTildeVec = r2TildeVec;
		rVec = r2Vec;
	end	

	if its > maxIT
		warning('Exceed the maximum iterate numbers');
		disp(['The iterative process stops at residual = ' sprintf('%10.4f',resnorm)]);		
	end
	if 2==nargout, varargout{1} = its; end
	clear rVec rTildeVec pVec p2Vec
% tMtV_
% tPtV_	
end

