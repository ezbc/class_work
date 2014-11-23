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
        
scatter(A(1,:),A(2,:),'*')
axis([-1 1 -1 1])
axis('square')
set(gca,'Fontsize',16)
set(gca,'Xtick',[-1 0 1])
set(gca,'Ytick',[-1 0 1])
print -djpeg cluster_data.jpg
