% spectral clustering
clear
close all

r=2; % spectral emdedding dimension
k=3; % k in kmeans

% generate synthetic sim matrix with clusters
k=4;
m(1:k) = ceil(10*rand(k,1));
n= sum(m);
P = rand(n,n)/3;
for i=1:k
    off = sum(m(1:i))-m(i);
    for j=1:m(i)
        for L=j:m(i)
            ii = j+off;
            jj = L+off;
            P(ii,jj) = 1 + rand/3;
            P(jj,ii) = P(ii,jj);
            if ii==jj, P(ii,jj)=1; end
        end
    end
end


% shuffled sim matrix to simulate unknown clustering
n = length(P);
[p,q]= sort(randn(n,1));
S = P(q,q);

% weight matrix from sims
W = exp(S); W = W/max(max(W)); 

% degree matrix
for i=1:n
    D(i,i) = sum(W(i,:));
end
% laplacian
L = D-W;

% eigendecomp to get embedding
[T,E,V] = svd(L);
U = T(:,n-r:n-1);  

% kmeans++ to cluster rows of U
[L_C,C] = kmeans(U',k);

% sort similarities according to clustering
[p,kq] = sort(L_C);
kS = S(kq,kq);

% display results
figure(1)
subplot(1,3,1)
imagesc(P)
axis('square')
title('ideal sim matrix')
subplot(1,3,2)
imagesc(W)
axis('square')
title('rand scrambled sims')
subplot(1,3,3)
imagesc(kS)
axis('square')
title('spectral kmeans')



