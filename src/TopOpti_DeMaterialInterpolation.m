function DeE = TopOpti_DeMaterialInterpolation(xPhys)
	global modulusMin_;
	global penalty_;	
	global modulus_;
	DeE = penalty_*(modulus_-modulusMin_)' .* xPhys.^(penalty_-1);
end