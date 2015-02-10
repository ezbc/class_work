feature(i).b=zeros(length(feature(i).A(:,1)),1)-1;

for j=[3]
    
    feature(i).b(j,1)=1;
    
end

find(feature(i).b==1)
save('feature_hough.mat','feature');