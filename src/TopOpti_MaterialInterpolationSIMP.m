function E = TopOpti_MaterialInterpolationSIMP(xPhys)
	global modulusMin_;
	global SIMPpenalty_;	
	global modulus_;
	E = modulusMin_+double(xPhys(:)').^SIMPpenalty_ .* (modulus_-modulusMin_);
end