function varargout = Solving_CG_GMGS(printP)
	global U_; 
	global F_; 
	global tol_; 
	global maxIT_; 
	
	%U_ = Solving_PreconditionedConjugateGradientSolver_previous(@Solving_KbyU_MatrixFree_previous, @Solving_Vcycle_previous, F_, tol_, maxIT_, printP, U_);
	[U_, its] = Solving_PreconditionedConjugateGradientSolver(@Solving_KbyU_MatrixFree, @Solving_Vcycle, F_, tol_, maxIT_, printP, U_);
	if 1==nargout, varargout{1} = its; end
end
