function xPhys = TopOpti_DualizeDesignVariable(xTilde)
	global beta_;
	global eta_;
	xPhys = (tanh(beta_*eta_) + tanh(beta_*(xTilde-eta_))) / (tanh(beta_*eta_) + tanh(beta_*(1-eta_)));	
end