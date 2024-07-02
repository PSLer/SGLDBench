function E = TopOpti_MaterialInterpolationSIMP(xPhys)
	global modulusMin_;
	global SIMPpenalty_;	
	global modulus_;
	E = modulusMin_+xPhys(:)'.^SIMPpenalty_ .* (modulus_-modulusMin_);
end