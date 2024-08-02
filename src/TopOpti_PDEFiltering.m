function tar = TopOpti_PDEFiltering(src)
	global TF_;
	global LF_;	
	global Permut_;		
	global KF_;  
	global maxIT_;
    global PDEfilterSolver_;
    src = double(src);
	AtV = @(x) KF_*x;
	PtV = @(x) LF_'\(LF_\x);				
	tar = TF_' * Solving_CG_standard(AtV, PtV, TF_*src, 1.0e-6, maxIT_, 'printP_OFF');	
end