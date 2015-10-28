% clear all;
close all;

load feature_hough.mat
upval=10;


A=[];
b=[];


count=1;
for i=1:length(feature)
    for j=1:length(feature(i).H_subimage)
    
        hough_feature1.hough_max_SHG(count,1)=max(max(feature(i).H_subimage(j).subimage));
        hough_feature1.hough_avg_SHG(count,1)=mean(mean(feature(i).H_subimage(j).subimage));
        hough_feature1.hough_nonzero_avg_SHG(count,1)=mean(feature(i).H_subimage(j).subimage(find(feature(i).H_subimage(j).subimage~=0)));
        count=count+1;
    end
     b=[b; feature(i).b];
end

A=[hough_feature1.hough_max_SHG hough_feature1.hough_avg_SHG hough_feature1.hough_nonzero_avg_SHG];

count =1;
for i=1:9:length(feature)-6
   b_temp=b(i:i+8);
   if(sum(b_temp)==-9)
       b_new(count,1)=-1;
   else b_new(count,1)=1;
   end
   count=count+1;
end

b=b_new;
