function tar = TopOpti_PDEFiltering(src)
	global TF_;
	global LF_;	
	global KF_;  
	global maxIT_;
    src = double(src);
	AtV = @(x) KF_*x;
	PtV = @(x) LF_'\(LF_\x);
	src = TF_*src;
	tar = Solving_CG_standard(AtV, PtV, src, 1.0e-6, maxIT_, 'printP_OFF');
	tar = TF_' * tar;
end