% create data for clustering
clear
close all

n=100;

cnt=1;
while cnt < n
    a = 2*rand(2,1)-1;
    if norm(a) < 1/3
        A(:,cnt) = a;
        cnt=cnt+1;
    elseif (norm(a) > .5 & norm(a) < .6)
        A(:,cnt) = a;
        cnt=cnt+1;
    end
end
        
figure(1);

[L,C] = kmeans(A,2);
hold on;
scatter(A(1,:),A(2,:),'*')
scatter(A(1,L == 1),A(2,L == 1),'go')
scatter(A(1,L == 2),A(2,L == 2),'bo')
axis([-1 1 -1 1])
axis('square')
set(gca,'Fontsize',16)
set(gca,'Xtick',[-1 0 1])
set(gca,'Ytick',[-1 0 1])
hold off;

r=2; % spectral emdedding dimension
%construct similarity matrix
m = length(A);
S = zeros(m,m);
for i=1:m
    for j=1:m
        if i==j
            S(i,j) = 0;
        else
            S(i,j) = exp(-norm(A(:,i)-A(:,j))^2/2) ;
        end
    end
end

S = S/max(max(S)); 
% degree matrix
for i=1:m
    D(i,i) = sum(S(i,:));
end
% laplacian
L = sqrt(D)\S/sqrt(D);
%L = D-S;
% eigendecomp to get embedding
[T,E,V] = svd(L);
U = T(:,m-r:m-1);  
% kmeans++ to cluster rows of U
[L,C] = kmeans(U',2);
[p,kq] = sort(L);
kS = S(kq,kq);

figure(2);

hold on;
scatter(A(1,:),A(2,:),'*')
scatter(A(1,L == 1),A(2,L == 1),'go')
scatter(A(1,L == 2),A(2,L == 2),'bo')
axis([-1 1 -1 1])
axis('square')
set(gca,'Fontsize',16)
set(gca,'Xtick',[-1 0 1])
set(gca,'Ytick',[-1 0 1])
hold off;

