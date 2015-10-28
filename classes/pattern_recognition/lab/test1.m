% clear all;
close all;

%load feature_hough





i=19;



%^%%%%%%%%%%%%%%%%%%%%%% CHANGE I%%%%%%%%%%%%%%%%%%%%%%
HE_IMAGE='TransformedTumour Pancreas 84-1.tif';
SHG_IMAGE='pancreatic tumor 2000 - 84-1.tif';
% HE_IMAGE='G3 #145 tumor pancreas.tif';
% SHG_IMAGE='SHGresized.tif';



[A,H_subimage,H_block]=hough_feature_generator(HE_IMAGE,SHG_IMAGE,50);

figure
imshow(H_subimage)

% test_feature=hough_feature_generator(HE_IMAGE,SHG_IMAGE,50);


