function vFpassiveElements = TopOpti_GetVolumeFractionOfPassiveElements()
	global meshHierarchy_;
	global passiveElements_; 
	vFpassiveElements = numel(passiveElements_) / meshHierarchy_(1).numElements;
end