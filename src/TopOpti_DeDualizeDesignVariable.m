function dx = TopOpti_DeDualizeDesignVariable(xTilde)
	global beta_;
	global eta_;
	dx = beta_*(1-tanh(beta_*(xTilde-eta_)).*tanh(beta_*(xTilde-eta_)))/(tanh(beta_*eta_)+tanh(beta_*(1-eta_)));		
end