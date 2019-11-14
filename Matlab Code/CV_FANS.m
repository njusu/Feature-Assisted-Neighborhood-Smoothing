function [L2, L1, L22, L11] = CV_FANS(A, X, LAMBDA, ITER)
nlambda = length(LAMBDA);
PE1 = zeros(ITER, nlambda); PE2 = zeros(ITER, nlambda);
PE3 = zeros(ITER, nlambda); PE4 = zeros(ITER, nlambda);
n = size(A,1);
% calculate distance s, don't need d0 since d0_train will be calculated
% after splitting data to training and validation set
s = pdist2(X,X);

parfor iter=1:ITER
    ind_valid = sort( datasample(1:n, round(n/10), 'Replace', false) ); % or use  randperm(n,round(n/10))
    % can sort or not sort, order doesn't matter
    ind_train = 1:n; ind_train(ind_valid) = [];
    A_train = A(ind_train,ind_train);
    X_train = X(ind_train,:);
    X_valid = X(ind_valid,:);
    
    % Fit training model
    n_train = length(ind_train);
    d0_train = zeros(n_train, n_train);
	A_sq_train = A_train*A_train;
	for i = 1:(n_train-1)
		for j = (i+1):n_train
			temp = abs(A_sq_train(i, :) - A_sq_train(j, :));
            temp([i,j])=[];
            d0_train(i, j) = ( max(temp) +rand(1,1)/n_train )/n_train;  
			d0_train(j, i) = d0_train(i, j);
		end
    end
    s_train = s(ind_train,ind_train);
    
    for k=1:nlambda
        lambda = LAMBDA(k);
        d_train = d0_train./max(d0_train(:)) + lambda * s_train./max(s_train(:));
        % Define neighborhood 
        Kernel_mat = zeros(n_train, n_train);
        h_train = sqrt(log(n_train)/n_train);
        for i = 1:n_train
            di = d_train(i,:); di(i)=[];
		    Kernel_mat(i, :) =  (d_train(i, :) <= quantile(di, h_train));
            if(sum(Kernel_mat(i, :)) > 1)
                Kernel_mat(i, i) = 0;
            end
        end
        
        Kernel_mat = Kernel_mat ./ (sum(Kernel_mat, 2)*ones(1, n_train));  % Each row has been normalized(under L1)
	    W_train = Kernel_mat * A_train;
	    P_hat_train = (W_train + W_train') / 2;
        
        % Find the nearest nodes to the validation set
        index = zeros(1,length(ind_valid));
        for i = 1:length(ind_valid)
            [~,index(i)] = min(s(ind_valid(i),ind_train));
        end
        PE1(iter,k) = mean(mean( (P_hat_train(index,:)-A(ind_valid,ind_train)).^2 ));
        PE2(iter,k) = mean(mean( abs(P_hat_train(index,:)-A(ind_valid,ind_train)) ));
        temp = abs(P_hat_train(index,:)-A(ind_valid,ind_train)); temp = sort(temp(:));
        PE3(iter,k) = mean( temp(1:(0.09*n^2/4)).^2 );
        PE4(iter,k) = mean( temp(1:(0.09*n^2/4)) );
    end
    fprintf('iter = %3g \n', iter);
end
L2 = mean(PE1,1); L1 = mean(PE2,1);
L22 = mean(PE3,1); L11 = mean(PE4,1);
end

    
    
