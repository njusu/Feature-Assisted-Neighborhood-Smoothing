function [W_hat] = NBS(A,c)
	% INPUT:
	% A: adjacency matrix
	% OUTPUT:
	% W_hat: estimated probaility matrix
	% COPYRIGHT:
	% Estimating network edge probabilities by neighborhood smoothing
	% Yuan Zhang, Elizaveta Levina and Ji Zhu
	% Contact: Yuan Zhang yzhanghf@gmail.com
	
	N = size(A, 1);
	h = c*sqrt(log(N)/N);
	
	% Compute disimilarity measures
	D = zeros(N, N);
	A_sq = A*A;
	for i = 1:(N-1)
		for j = (i+1):N
			temp = abs(A_sq(i, :) - A_sq(j, :));
            temp([i,j])=[];
            D(i, j) = ( max(temp) +rand(1,1)/N )/N; 
			D(j, i) = D(i, j);
		end
	end
	
	for i = 1:N
        Di = D(i,:); Di(i)=[];
		Kernel_mat(i, :) =  (D(i, :) <= quantile(Di, h));
        if(sum(Kernel_mat(i, :)) > 1)
            Kernel_mat(i, i) = 0;
        end
	end
	
	Kernel_mat = Kernel_mat ./ (sum(Kernel_mat, 2)*ones(1, N)); 
	
	W_hat = Kernel_mat * A;
	W_hat = (W_hat + W_hat') / 2;
	
end
