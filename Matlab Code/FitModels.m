A = textread('C:/Users/Yi/Desktop/Methods Comparison/n=500/A_g1.txt');
A_2 = textread('C:/Users/Yi/Desktop/Methods Comparison/n=500/A2_g2.txt');
P_GT = textread('C:/Users/Yi/Desktop/Methods Comparison/n=500/P_g2.txt'); % for SBA

A = textread('C:/Users/Yi/Desktop/A_comm10.txt');
A_2 = textread('C:/Users/Yi/Desktop/A_real.txt');

% NBS
P_NBS =  NBS(A,1);

% SAS
[M0 N0] = size(A);
n = M0;
h = round(log(n));
d = mean(A);
[~, pos] = sort(d,'descend');
% estimate SAS
P_SAS(pos,pos) = SAS(A);

% SBA
G(:,:,1) = A;
G(:,:,2) = A_2;
% run Demo_crossvalidation.m to choose Delta
Delta = 0.1; % n=500, all g can use 0.1, g4 can also try 0.05 
B = estimate_blocks_directed(G,Delta);
[~,P_SBA] = histogram3D(G,B); % ~ was H

% USVT
P_USVT = Method_chatterjee((A-0.5).*2);
P_USVT = P_USVT/2+0.5;

% write results out to R
dlmwrite('P_NBS.txt',P_NBS);
dlmwrite('P_SAS.txt',P_SAS);
dlmwrite('P_SBA.txt',P_SBA);
dlmwrite('P_USVT.txt',P_USVT);
dlmwrite('P_NBSF.txt',P_NBSF);
dlmwrite('P_true.txt',P);

%colormap jet(1000);
imagesc(P_NBS);
imagesc(P_SAS);
imagesc(P_SBA);
imagesc(P_USVT);