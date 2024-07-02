function D = FEA_ElementElasticityMatrix(S)
	D = zeros(48);
	for ii=1:8
		index = (ii-1)*6+1:ii*6;
		D(index,index) = S;
	end
end
