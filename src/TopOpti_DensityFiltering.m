function tar = TopOpti_DensityFiltering(src, opt)
	global H_;
	global Hs_;	
	if 1==opt
		tar = H_*(src./Hs_);
	elseif 0==opt
		tar = H_*src./Hs_;
	else
		error('Wrong option for checker board filtering')
	end
end