% HW 9 ECE 532

% ------------------------------------------------------------------------------
% create data for clustering
% ------------------------------------------------------------------------------
clear
close all

n=100;
sigma=2;
r=1;
k=2;
cluster_sep_scale = 3; % higher means cluster are more separated

cnt=1;
while cnt < n
    a = sigma*rand(2,1)-1;
    if norm(a) < 1/cluster_sep_scale
        A(:,cnt) = a;
        cnt=cnt+1;
    elseif (norm(a) > .5 & norm(a) < .6)
        A(:,cnt) = a;
        cnt=cnt+1;
    end
end

% ------------------------------------------------------------------------------
% Cluster with kmeans only
% ------------------------------------------------------------------------------

[L_k,C] = kmeans(A,2);

figure(1);
hold on;
gscatter(A(1,:),A(2,:),L_k)
axis([-1 1 -1 1])
title('Regular Kmeans')
axis('square')
set(gca,'Fontsize',16)
set(gca,'Xtick',[-1 0 1])
set(gca,'Ytick',[-1 0 1])
hold off;

% ------------------------------------------------------------------------------
% Cluster spectrally
% ------------------------------------------------------------------------------
P = A';
n = length(P);

S = P;

% Create weight matrix
W = zeros(length(P), length(P));
for i=1:length(P)
    for j=1:length(P)
        %W(i,j) = pdist2(P(i,:), P(j,:));
        if i==j
            W(i,j) = 0;
        else
            %W(i,j) = exp(-norm(P(i,:) - P(j,:))^2/2);
            %W(i,j) = exp(-norm(P(i,:) - P(j,:)));
            %W(i,j) = exp(-sum((P(i,:).^2 - P(j,:).^2).^0.5)/(2*sigma^2));
            ri = (P(i,1)^2 + P(i, 2)^2)^0.5;
            rj = (P(j,1)^2 + P(j, 2)^2)^0.5;
            W(i, j) = exp(-norm(ri - rj)^2/2);
        end
    end
end
%W = exp(W); 
W = W/max(max(W)); % normalize weight matrix

% degree matrix
for i=1:n
    D(i,i) = sum(W(i,:));
end

% laplacian
%L = D-W;
L = sqrt(D)\W/sqrt(D);

% eigendecomp to get embedding
[T,E,V] = svd(L);
U = T(:,n-r:n-1);

% kmeans++ to cluster rows of U
[L_C,C] = kmeans(U',k);

% sort similarities according to clustering
[p,kq] = sort(L_C);
kS = W(kq,kq);

%disp('kS')
%disp(kS)

% display results
figure(2);
% Original weight matrix
subplot(1,3,1)
imagesc(W)
axis('square')
title('rand scrambled sims')

% Clustered matrix
subplot(1,3,2)
imagesc(kS)
axis('square')
title('spectral kmeans')

% Clustered dataspace
subplot(1,3,3)
%gscatter(P(kq,1),P(kq,2),p)
gscatter(P(:,1),P(:,2),L_C)
%scatter(P(L_C==1,1),P(L_C==1,2),'go')
%scatter(P(L_C==2,1),P(L_C==2,2),'bo')
axis([-1 1 -1 1])
axis('square')
set(gca,'Fontsize',16)
set(gca,'Xtick',[-1 0 1])
set(gca,'Ytick',[-1 0 1])
print -djpeg cluster_data.jpg




