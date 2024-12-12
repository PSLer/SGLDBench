function S = FEA_HookeLaw()
	global modulus_; 
	global poissonRatio_;
	E = 1; nu = poissonRatio_;
	S = zeros(6);
	cons1 = (1+nu)*(1-2*nu);
	cons2 = 2*(1+nu);
	S(1,1) = E*(1-nu)/cons1; S(1,2) = E*nu/cons1; S(1,3) = S(1,2);
	S(2,1) = S(1,2); S(2,2) = S(1,1); S(2,3) = S(2,1);
	S(3,1) = S(1,3); S(3,2) = S(2,3); S(3,3) = S(2,2);
	S(4,4) = E/cons2;
	S(5,5) = S(4,4);
	S(6,6) = S(5,5);
end