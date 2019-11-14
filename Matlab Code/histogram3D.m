function [H,P] = histogram3D(G,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function histogram3D(G,B)
% computes the 3D histogram H and the corresponding Probability matrix P
% from the cluters B
%
% Input: G, n x n x T graph
%        B, clusters
% Output: H, 3D histogram
%         P, the corresponding probability matrix of size nxn
%
% Stanley Chan @ Harvard
% Feb 13, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get n
n = size(G,1);
T = size(G,3);

% Initialization
P = zeros(n,n);
H = zeros(length(B), length(B));

% Loop through all the clusters
for ki=1:length(B)
    for kj=1:length(B)
        % Obtain the indices in cluster ki and kj
        I = B{ki};
        J = B{kj};
        
        % Compute the 3D histogram
        H(ki,kj) = sum(sum(sum(G(I,J,:))))/(T*length(I)*length(J));
        
        % Loop through the indices in cluster I and J
        % Compute the corresponding probability matrix
        for i1=1:length(I)
            for j1=1:length(J)
                vi = I(i1);
                vj = J(j1);
                P(vi,vj) = H(ki,kj);
            end
        end     
    end
end
