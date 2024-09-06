#ifndef __MMAseq__
#define __MMAseq__
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

class MMAseq{

public:
	MMAseq(int n,int m);
	MMAseq(int n,int m, double a, double c, double d);
	~MMAseq();

	void SetAsymptotes(double init, double decrease, double increase);
	
	void ConstraintModification(bool conMod);
	
	void Update(double *xval, double *dfdx, double *gx, double *dgdx, double *xmin, double *xmax);

	void Reset(){ k=0; };
	double *xo1,*xo2;
private:

	int n,m,k;

	double asyminit,asymdec,asyminc;

	bool constraintModification; 
  	
	double *a, *c, *d;

	double *y,z;

	double *lam,*mu,*s;

	double *L,*U,*alpha,*beta,*p0,*q0,*pij,*qij,*b,*grad,*Hess;

	//double *xo1,*xo2;

	void GenSub(double *xval, double *dfdx, double *gx, double *dgdx, double *xmin, double *xmax);

	void SolveDSA(double *x);
	void SolveDIP(double *x);


	void XYZofLAMBDA(double *x);

	void DualGrad(double *x);
	void DualHess(double *x);
	void DualLineSearch();
	double DualResidual(double *x, double epsi);

	void Factorize(double *K, int n);
	void Solve(double *K, double *x, int n);
	// Math helpers
	double Min(double d1, double d2);
	double Max(double d1, double d2);
	int Min(int d1, int d2);
	int Max(int d1, int d2);
	double Abs(double d1);
};

#endif
