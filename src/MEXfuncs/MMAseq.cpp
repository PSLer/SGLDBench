#include "MMAseq.h"
#include <iostream>
#include <math.h>

/*

BETA VERSION  0.99

Sequential MMA solver using a dual interior point method

Code by Niels Aage, February 2013

The class solves a general non-linear programming problem 
on standard from, i.e. non-linear objective f, m non-linear 
inequality constraints g and box constraints on the n
design variables xmin, xmax.

       min_x^n f(x)
       s.t. g_j(x) < 0,   j = 1,m
       xmin < x_i < xmax, i = 1,n

Each call to Update() sets up and solve the following 
convex subproblem:

  min_x     sum(p0j./(U-x)+q0j./(x-L)) + a0*z + sum(c.*y + 0.5*d.*y.^2)
 
  s.t.      sum(pij./(U-x)+qij./(x-L)) - ai*z - yi <= bi, i = 1,m
            Lj < alphaj <=  xj <= betaj < Uj,  j = 1,n
            yi >= 0, i = 1,m        
            z >= 0.
 
NOTE: a0 == 1 in this implementation !!!!

*/

// PUBLIC
MMAseq::MMAseq(int nn,int mm){
	n = nn;
	m = mm;
	
	asyminit = 0.5;//0.2;
	asymdec = 0.7;//0.65;
	asyminc = 1.2;//1.08;
	
	constraintModification = false;
	
	k = 0;

	a = new double[m];
	c = new double[m];
	d = new double[m];

	for (int i=0;i<m;i++){
		a[i] = 0.0;
		c[i] = 1000.0;
		d[i] = 0.0;
	}

	y = new double[m];
	lam = new double[m];

	L = new double[n];
	U = new double[n];

	alpha = new double[n];
	beta = new double[n];
	
	p0 = new double[n];
	q0 = new double[n];
	pij = new double[n*m];
	qij = new double[n*m];
	b = new double[m];

	xo1 = new double[n];
	xo2 = new double[n];

	grad = new double[m];
	mu = new double[m];
	s = new double[2*m];

	Hess = new double[m*m];
}

MMAseq::MMAseq(int nn,int mm, double ai, double ci, double di){
	n = nn;
	m = mm;
	
	asyminit = 0.5;//0.2;
	asymdec = 0.7;//0.65;
	asyminc = 1.2;//1.08;
	
	constraintModification = false;
	
	k = 0;

	a = new double[m];
	c = new double[m];
	d = new double[m];

	for (int i=0;i<m;i++){
		a[i] = ai;
		c[i] = ci;
		d[i] = di;
	}

	y = new double[m];
	lam = new double[m];

	L = new double[n];
	U = new double[n];

	alpha = new double[n];
	beta = new double[n];
	
	p0 = new double[n];
	q0 = new double[n];
	pij = new double[n*m];
	qij = new double[n*m];
	b = new double[m];

	xo1 = new double[n];
	xo2 = new double[n];

	grad = new double[m];
	mu = new double[m];
	s = new double[2*m];

	Hess = new double[m*m];


}

MMAseq::~MMAseq(){
	delete [] a;
	delete [] b;
	delete [] c;
	delete [] d;
	delete [] y;
	delete [] lam;
	delete [] L;
	delete [] U;
	delete [] alpha;
	delete [] beta;
	delete [] p0;
	delete [] q0;
	delete [] pij;
	delete [] qij;
	delete [] xo1;
	delete [] xo2;
	delete [] grad;
	delete [] mu;
	delete [] s;
	delete [] Hess;
}

void MMAseq::SetAsymptotes(double init, double decrease, double increase){
	
	// asymptotes initialization and increase/decrease
	asyminit = init;
	asymdec = decrease;
	asyminc = increase;
}

void MMAseq::ConstraintModification(bool conMod){
	
	constraintModification = conMod;
  
}

void MMAseq::Update(double *xval, double *dfdx, double *gx, double *dgdx, double *xmin, double *xmax){

	// Generate the subproblem
	GenSub(xval,dfdx,gx,dgdx,xmin,xmax);
	
	// Update xolds
#pragma omp parallel for
    for (int i=0;i<n;i++){
        xo2[i]=xo1[i];
        xo1[i]=xval[i];
    }

	// Solve the dual with an interior point method
	 SolveDIP(xval);

	// Solve the dual with a steepest ascent method
	//SolveDSA(xval);

}

// PRIVATE
void MMAseq::SolveDIP(double *x){

	for (int j=0;j<m;j++){
                lam[j]=c[j]/2.0;
		mu[j]=1.0;
        }

	double tol = 1.0e-9*sqrt(m+n);
	double epsi=1.0;
        double err = 1.0;
        int loop;

	while (epsi>tol){
		
		loop = 0;
		while (err>0.9*epsi && loop<100){
			loop++;

			// Set up newton system
			XYZofLAMBDA(x);
			DualGrad(x);
			for (int j=0;j<m;j++){	
				grad[j] = -1.0*grad[j] -epsi/lam[j];
			}
			DualHess(x);
			// Solve Newton system
			if (m>1){ 
				Factorize(Hess,m);
				Solve(Hess,grad,m);
				for (int j=0;j<m;j++){
                 	s[j]=grad[j];
           	}
			}
			else{
				s[0] = grad[0]/Hess[0];
			}

			// Get the full search direction
			for (int i=0;i<m;i++){
     		    s[m+i]= -mu[i]+epsi/lam[i]-s[i]*mu[i]/lam[i];
       	}
		
			// Perform linesearch and update lam and mu
			DualLineSearch();

			XYZofLAMBDA(x);

			// Compute KKT res
			err = DualResidual(x,epsi);			

		}
		epsi=epsi*0.1;
	}

}


void MMAseq::SolveDSA(double *x){

	for (int j=0;j<m;j++){
		lam[j]=1.0;
	}

	double tol = 1.0e-9*sqrt(m+n);
	double err = 1.0;
	int loop = 0;	

	while (err>tol && loop<500){
		loop++;
		XYZofLAMBDA(x);
		DualGrad(x);	
		double theta = 1.0;
		err=0.0;
		for (int j=0;j<m;j++){
			lam[j] = Max(0.0,lam[j]+theta*grad[j]);
			err+=grad[j]*grad[j];
		}	
		err=sqrt(err);
	}
}

double MMAseq::DualResidual(double *x, double epsi){

	double *res = new double[2*m];
	
        for (int j=0;j<m;j++){
                res[j]= - b[j] - a[j]*z - y[j]+mu[j];
		res[j+m]= mu[j]*lam[j] - epsi;
                for (int i=0;i<n;i++){
                        res[j]+=pij[i*m+j]/(U[i]-x[i]) + qij[i*m+j]/(x[i]-L[i]);
                }
        }

	double nrI=0.0;
	for (int i=0;i<2*m;i++){
		if (nrI<Abs(res[i])){
			nrI=Abs(res[i]);
		}
	}

	delete [] res;

	return nrI;
}

void MMAseq::DualLineSearch(){
	
	double theta=1.005;
	for (int i=0;i<m;i++){
		if (theta < -1.01*s[i]/lam[i]){
			theta = -1.01*s[i]/lam[i];
		}
		if (theta < -1.01*s[i+m]/mu[i]){
                        theta = -1.01*s[i+m]/mu[i];
                }
	}
	theta = 1.0/theta;
	
	for (int i=0;i<m;i++){
		lam[i] = lam[i] + theta*s[i];
		mu[i] = mu[i] + theta*s[i+m];
	}
	
}

void MMAseq::DualHess(double *x){

	double *df2 = new double[n];
	double *PQ = new double[n*m];
#pragma omp parallel for	
	for (int i=0;i<n;i++){
		double pjlam, qjlam;
		pjlam = p0[i];
        qjlam = q0[i];
        for (int j=0;j<m;j++){
            pjlam+=pij[i*m+j]*lam[j];
            qjlam+=qij[i*m+j]*lam[j];
			PQ[i*m+j] = pij[i*m+j]/pow(U[i]-x[i],2.0) - qij[i*m+j]/pow(x[i]-L[i],2.0);
        }
		df2[i] = -1.0/(2.0*pjlam/pow(U[i]-x[i],3.0) + 2.0*qjlam/pow(x[i]-L[i],3.0));
		double xp = (sqrt(pjlam)*L[i] + sqrt(qjlam)*U[i]) / (sqrt(pjlam)+sqrt(qjlam));
		if (xp<alpha[i]){
			df2[i]=0.0;
		}
		if (xp>beta[i]){
			df2[i]=0.0;
		}
	}
	
	// Create the matrix/matrix/matrix product: PQ^T * diag(df2) * PQ
	double *tmp = new double[n*m];
	for (int j=0;j<m;j++){
#pragma omp parallel for		
		for (int i=0;i<n;i++){
			tmp[j*n+i]=0.0;
			tmp[j*n+i] += PQ[i*m+j]*df2[i];
			/*
			for (int k=0;k<n;k++){
				if (k==i){
					tmp[j*n+i] += PQ[k*m+j]*df2[k];
				}
			}
			*/

		}
	}

	for (int i=0;i<m;i++){
                for (int j=0;j<m;j++){
                        Hess[i*m+j]=0.0;
                        for (int k=0;k<n;k++){
				Hess[i*m+j]+=tmp[i*n+k]*PQ[k*m+j];	
			}
		}
	}

	double lamai=0.0;
	for (int j=0;j<m;j++){
		if (lam[j]<0.0){
			lam[j]=0.0;
		}
		lamai+=lam[j]*a[j];
		if (lam[j]>c[j]){
			Hess[j*m+j]+= -1.0;
		}
		Hess[j*m+j]+= -mu[j]/lam[j]; 
	}
	
	if (lamai>0.0){
		for (int j=0;j<m;j++){
			for (int k=0;k<m;k++){
				Hess[j*m+k]+= -10.0*a[j]*a[k];
			}
		}
	}

	// pos def check
	double HessTrace=0.0;
	for (int i=0;i<m;i++){
        	HessTrace+=Hess[i*m+i];
	}
    	double HessCorr = 1e-4*HessTrace/m;

	if (-1.0*HessCorr<1.0e-7){
		HessCorr= -1.0e-7;
	}

	for (int i=0;i<m;i++){
                Hess[i*m+i]+=HessCorr;
        }


	delete [] df2;
	delete [] PQ;
	delete [] tmp;

}

void MMAseq::DualGrad(double *x){
	for (int j=0;j<m;j++){
		grad[j]= - b[j] - a[j]*z - y[j];
		for (int i=0;i<n;i++){
			grad[j]+=pij[i*m+j]/(U[i]-x[i]) + qij[i*m+j]/(x[i]-L[i]);
		}
	}
}

void MMAseq::XYZofLAMBDA(double *x){

	double lamai=0.0;
	for (int i=0;i<m;i++){
		if (lam[i]<0.0){
			lam[i]=0;
		}
		y[i] = Max(0.0,lam[i]-c[i]); //Note y=(lam-c)/d - however d is fixed at one !!
		lamai+=lam[i]*a[i];
	}
	z=Max(0.0,10.0*(lamai-1.0)); // SINCE a0 = 1.0

#pragma omp parallel for
	for (int i=0;i<n;i++){
		double pjlam, qjlam;
		pjlam = p0[i];
		qjlam = q0[i];
		for (int j=0;j<m;j++){
			pjlam+=pij[i*m+j]*lam[j];
			qjlam+=qij[i*m+j]*lam[j];
		}
		x[i] = (sqrt(pjlam)*L[i] + sqrt(qjlam)*U[i]) / (sqrt(pjlam)+sqrt(qjlam));
		if (x[i]<alpha[i]){
			x[i]=alpha[i];
		}
		if (x[i]>beta[i]){
			x[i]=beta[i];
		}

	}
	

}

void MMAseq::GenSub(double *xval, double *dfdx, double *gx, double *dgdx, double *xmin, double *xmax){

	

	// forward the iterator
	k++;

	// Set asymptotes
	if (k<3){
#pragma omp parallel for		
		for(int i=0;i<n;i++){
			L[i] = xval[i] - asyminit*(xmax[i] - xmin[i]);
			U[i] = xval[i] + asyminit*(xmax[i] - xmin[i]);
		}
	}
	else{
#pragma omp parallel for		
         for(int i=0;i<n;i++){
			double gamma, helpvar;
			helpvar = (xval[i]-xo1[i])*(xo1[i]-xo2[i]);
			if (helpvar < 0.0){
				gamma = asymdec;
			}
			else if (helpvar > 0.0){
				gamma = asyminc;
			}
			else{ gamma = 1.0;}			 
			L[i] = xval[i] - gamma*(xo1[i] - L[i]);
                        U[i] = xval[i] + gamma*(U[i] - xo1[i]);
			
			double xmi,xma;
			xmi = Max(1.0e-5,xmax[i]-xmin[i]);			
			L[i]=Max(L[i],xval[i]-100.0*xmi);
			L[i]=Min(L[i],xval[i]-1.0e-5*xmi);
			U[i]=Max(U[i],xval[i]+1.0e-5*xmi);
			U[i]=Min(U[i],xval[i]+100*xmi);

			xmi = xmin[i]-1.0e-6;
			xma = xmax[i]+1.0e-6;
			if (xval[i]<xmi)
			{
				L[i]=xval[i]-(xma-xval[i])/0.9;
				U[i]=xval[i]+(xma-xval[i])/0.9;
			}
			if (xval[i]>xma)
			{
				L[i]=xval[i]-(xval[i]-xmi)/0.9;
				U[i]=xval[i]+(xval[i]-xmi)/0.9;
			}

			
        }
	}

	// Set bounds and the coefficients for the approximation
	
	double feps = 1.0e-6;
#pragma omp parallel for	
	for(int i=0;i<n;i++){
		double dfdxp, dfdxm;
		alpha[i] = Max(xmin[i],0.9*L[i]+0.1*xval[i]);
		beta[i]  = Min(xmax[i],0.9*U[i]+0.1*xval[i]);

		dfdxp = Max(0.0,dfdx[i]);
		dfdxm = Max(0.0,-1.0*dfdx[i]);		
		
		// C++
		p0[i] = pow(U[i]-xval[i],2.0)*( dfdxp + 0.001*Abs(dfdx[i]) + 0.5*feps/(U[i]-L[i]));
                q0[i] = pow(xval[i]-L[i],2.0)*( dfdxm + 0.001*Abs(dfdx[i]) + 0.5*feps/(U[i]-L[i]));


		for (int j=0;j<m;j++){

			dfdxp = Max(0.0,dgdx[i*m+j]);
	                dfdxm = Max(0.0,-1.0*dgdx[i*m+j]);
			
			// ROBUST WAY FOR GENERAL CONSTRAINTS
			if (constraintModification){
				pij[i*m+j] = pow(U[i]-xval[i],2.0)*( dfdxp + 0.001*Abs(dgdx[i]) + 0.5*feps/(U[i]-L[i]));
				qij[i*m+j] = pow(xval[i]-L[i],2.0)*( dfdxm + 0.001*Abs(dgdx[i]) + 0.5*feps/(U[i]-L[i]));
			}else{
			// Fast for linear constraints
			 pij[i*m+j] = pow(U[i]-xval[i],2.0)*( dfdxp );
			 qij[i*m+j] = pow(xval[i]-L[i],2.0)*( dfdxm );
			}
		}
		

	}
	
	// The constant for the constraints
	for (int j=0;j<m;j++){
		b[j] = -gx[j];
		for (int i=0;i<n;i++){
			b[j]+= pij[i*m+j]/(U[i]-xval[i]) + qij[i*m+j]/(xval[i]-L[i]);
		}
	}	

}


void MMAseq::Factorize(double *K, int n){

	for (int s=0;s<n-1;s++){
		for (int i=s+1;i<n;i++){
			K[i*n+s] = K[i*n+s] / K[s*n+s];
			for (int j=s+1;j<n;j++){
				K[i*n+j] = K[i*n+j] - K[i*n+s]*K[s*n+j]; 				
			}
		}

	}	

}

void MMAseq::Solve(double *K, double *x, int n){

	for (int i=1;i<n;i++){
		double a = 0.0;
		for (int j=0;j<i;j++){
			a = a - K[i*n+j]*x[j];
		}
		x[i] = x[i] + a;
	}

	x[n-1] = x[n-1]/K[(n-1)*n+(n-1)];
	for (int i=n-2;i>=0;i--){
		double a=x[i];
		for (int j=i+1;j<n;j++){
			a = a - K[i*n+j]*x[j];
		}
		x[i] = a/K[i*n+i];
	}

}

double MMAseq::Min(double d1, double d2){
        return d1<d2 ? d1 : d2;
}

double MMAseq::Max(double d1, double d2){
        return d1>d2 ? d1 : d2;
}

int MMAseq::Min(int d1, int d2){
        return d1<d2 ? d1 : d2;
}

int MMAseq::Max(int d1, int d2){
        return d1>d2 ? d1 : d2;
}

double MMAseq::Abs(double d1){
        return d1>0 ? d1 : -1.0*d1;
}




