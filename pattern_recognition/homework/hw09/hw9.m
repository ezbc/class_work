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








