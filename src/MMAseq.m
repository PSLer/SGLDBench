%% This MMA implementation is a translation of the c-version developed by Niles Aage at the DTU
function [xnew, xold1, xold2] = MMAseq(mm, nn, xvalTmp, xmin, xmax, xold1, xold2, dfdx, gx, dgdx)
	global n; n = nn;
	global m; m = mm;
	
	global asyminit; asyminit = 0.5;%%0.2;
	global asymdec; asymdec = 0.7;%%0.65;
	global asyminc; asyminc = 1.2;%%1.08;	
	
	global k; k = 0;
	
	global a; a = zeros(m,1);
	global c; c = 100.0 * ones(m,1);
	global d; d = zeros(m,1);
	
	global y; y = zeros(m,1);
	global z; z = 0;
	global lam; lam = zeros(m,1);
	
	global L; L = zeros(1,n);
	global U; U = zeros(1,n);
	
	global alpha; alpha = zeros(1,n);
	global beta; beta = zeros(1,n);
	
	global p0; p0 = zeros(1,n);
	global q0; q0 = zeros(1,n);
	global pij; pij = zeros(m,n);
	global qij; qij = zeros(m,n); 
	global b; b = zeros(m,1);
	
	global xo1; xo1 = xold1;
	global xo2; xo2 = xold2;
	
	global grad; grad = zeros(m,1);
	global mu; mu = zeros(m,1);
	global s; s = zeros(2*m,1);
	global Hess; Hess = zeros(m,m);
	
	global xVal; xVal = xvalTmp(:)';
	
	Update(dfdx(:)', gx(:), reshape(dgdx,n,m)', xmin(:)', xmax(:)');
	xnew = xVal(:); xold1 = xo1(:); xold2 = xo2(:);
	
	clear -global n m asyminit asymdec asyminc k a c d y z lam L U alpha beta p0 q0 pij qij b xo1 xo2 grad mu s Hess xVal
end

function [xnew, xold1, xold2] = Update(dfdx, gx, dgdx, xmin, xmax)
	global xo1 xo2 xVal;
	
	%% Generate the subproblem
	GenSub(dfdx,gx,dgdx,xmin,xmax);	%%Checked
	%% Update xolds
	xo2 = xo1;
	xo1 = xVal;
	%% Solve the dual with an interior point method
	SolveDIP();
end

function GenSub(dfdx,gx,dgdx,xmin,xmax)
	global k asyminit L U alpha beta p0 q0 pij qij b xVal;
	
	%% forward the iterator
	k = 1;
	%% Set asymptotes
	L = xVal - asyminit*(xmax - xmin);
	U = xVal + asyminit*(xmax - xmin);
	%% Set bounds and the coefficients for the approximation
	feps = 1.0e-6;
	alpha = 0.9*L+0.1*xVal;
	tmpBool = alpha-xmin<0; alpha(tmpBool) = xmin(tmpBool);
	beta = 0.9*U+0.1*xVal;
	tmpBool = beta-xmax>0; beta(tmpBool) = xmax(tmpBool);
	
	dfdxp = dfdx; dfdxp(dfdxp<0.0) = 0.0;
	dfdxm = -1.0*dfdx; dfdxm(dfdxm<0.0) = 0.0;
	p0 = (U-xVal).^2.0 .*(dfdxp + 0.001*abs(dfdx) + 0.5*feps./(U-L));
	q0 = (xVal-L).^2.0 .*(dfdxm + 0.001*abs(dfdx) + 0.5*feps./(U-L));
	
	dfdxp = dgdx; dfdxp(dfdxp<0.0) = 0.0;
	dfdxm = -1.0*dgdx; dfdxm(dfdxm<0.0) = 0.0;
	pij = (U-xVal).^2.0 .* dfdxp;
	qij = (xVal-L).^2.0 .* dfdxm;	
	%% The constant for the constraints
	b = -gx + sum(pij./(U-xVal) + qij./(xVal-L), 2);
end

function SolveDIP()
	global n m lam mu c grad s Hess;
	
	lam = c/2.0;
	mu = ones(size(mu));
	tol = 1.0e-9*sqrt(m+n);
	epsi = 1.0;
	err = 1.0;
	while epsi > tol
		loop = 0;
		while err>0.9*epsi && loop<100
			loop = loop + 1;	
			%% Set up newton system
			XYZofLAMBDA();		
			DualGrad();		
			for jj=1:m
				grad(jj) = -1.0 * grad(jj) - epsi/lam(jj);
			end
			DualHess();
			%% Solve Newton system
			s(1:m,1) = Hess\grad;
			%% Get the full search direction
			s(m+1:2*m,1) = -mu + epsi./lam(:) - s(1:m,1).*mu(:)./lam(:);
			%% Perform linesearch and update lam and mu
			DualLineSearch();
			XYZofLAMBDA();	
			%% Compute KKT res
			err = DualResidual(epsi);			
		end
		epsi=epsi*0.1;
	end
end

function XYZofLAMBDA()
	global lam a c y z p0 q0 U L alpha beta pij qij xVal;

	lam(lam<0.0) = 0;
	y = lam - c; y(y<0.0) = 0.0;
	lamai = lam(:)' * a;	
	z = max(0.0, 10.0*(lamai-1.0));	
	
	pjlam = p0 + sum(pij.*lam, 1);
	qjlam = q0 + sum(qij.*lam, 1);
	xVal = (sqrt(pjlam).*L + sqrt(qjlam).*U) ./ (sqrt(pjlam)+sqrt(qjlam));
	tmpBool = xVal-alpha<0; xVal(tmpBool) = alpha(tmpBool);
	tmpBool = xVal-beta>0; xVal(tmpBool) = beta(tmpBool);
	
	clear pjlam qjlam
end

function grad = DualGrad()
	global a y U L grad b z pij qij xVal;
	
	grad = -b - a*z - y;
	grad = grad + sum(pij./(U-xVal) + qij./(xVal-L),2);
end

function DualHess()
	global m p0 q0 pij qij U L lam alpha beta Hess a c mu xVal;
	
	pjlam = p0 + sum(pij.*lam, 1);
	qjlam = q0 + sum(qij.*lam, 1);
	PQ = pij./(U-xVal).^2.0 - qij./(xVal-L).^2.0;
	df2 = -1.0 ./(2.0*pjlam./(U-xVal).^3.0 + 2.0*qjlam./(xVal-L).^3.0);
	xp = (sqrt(pjlam).*L + sqrt(qjlam).*U) ./ (sqrt(pjlam)+sqrt(qjlam));
	df2(xp-alpha<0) = 0.0;
	df2(xp-beta>0) = 0.0;
	%% Create the matrix/matrix/matrix product: PQ^T * diag(df2) * PQ
	tmp = PQ .* df2;
	for ii=1:m
		for jj=1:m
			Hess(jj,ii) = 0;
			Hess(jj,ii) = Hess(jj,ii) + tmp(ii,:)*PQ(jj,:)';
		end
	end
	lamai=0.0;
	for jj=1:m
		if lam(jj)<0.0, lam(jj) = 0.0; end
		lamai = lamai + lam(jj)*a(jj);
		if lam(jj)>c(jj)
			Hess(jj,jj) = Hess(jj,jj) -1.0;
		end
		Hess(jj,jj) = Hess(jj,jj) - mu(jj)/lam(jj); 
	end
	if lamai>0.0
		for jj=1:m
			for kk=1:m
				Hess(kk,jj) = Hess(kk,jj) - 10.0*a(jj)*a(kk);
			end
		end
	end	
	%% pos def check
	HessTrace = trace(Hess);
	HessCorr = 1e-4*HessTrace/m;
	if -1.0*HessCorr < 1.0e-7, HessCorr = -1.0e-7; end
	Hess = Hess + diag(repmat(HessCorr,m,1));
	
	clear pjlam qjlam df2 PQ tmp
end

function DualLineSearch()
	global m s lam mu;
	theta=1.005;
	for jj=1:m
		if theta < -1.01*s(jj)/lam(jj), theta = -1.01*s(jj)/lam(jj); end
		if theta < -1.01*s(jj+m)/mu(jj), theta = -1.01*s(jj+m)/mu(jj); end
	end
	theta = 1.0/theta;
	lam = lam + theta*s(1:m,:);
	mu = mu + theta*s(m+1:2*m,1);
end

function nrI = DualResidual(epsi)
	global m b a z lam mu U L pij qij y xVal;
	res = zeros(2*m,1);
	res(1:m,1) = -b - a.*z -y + mu;
	res(m+1:2*m,1) = mu.*lam - epsi;
	res(1:m,1) = res(1:m,1) + sum(pij./(U-xVal),2) + sum(qij./(xVal-L),2);
	
	nrI=0.0;
	for jj=1:2*m
		if nrI<abs(res(jj))
			nrI = abs(res(jj));
		end
	end
end