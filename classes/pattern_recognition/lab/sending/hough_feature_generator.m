%hough classificaition full code
% clear all;
% clc;
% close all

 function[A H_subimage H]=hough_feature_generator(HE_INPUT_IMAGE,SHG_INPUT_IMAGE,size_bwopen)

% SHG_INPUT_IMAGE='SHG1.tif';
% HE_INPUT_IMAGE='HE1.tif';




SHG_RESIZED_IMAGE='SHG_resized.tif';


%%%%main code started
a_mas=imread(SHG_INPUT_IMAGE);
a_HE=imread(HE_INPUT_IMAGE);



[m n f]=size(a_HE);
%
a_mas=imresize(a_mas,[m n]);

% imwrite(a_mas,SHG_RESIZED_IMAGE);

a_gray=rgb2gray(a_HE);

a_mas=im2double(a_mas);
a_SHG=a_mas;
bw=im2bw(a_mas,0.04);
bw=bwareaopen(bw,size_bwopen);

imshow(bw)

%%%creating intensity mask to remove white
a_gray=rgb2gray(a_HE);
temp=zeros([m n]);
for i=1:m
    for j=1:n
        if(a_gray(i,j)<190)
            temp(i,j)=a_gray(i,j);
        end
    end
end

temp= mat2gray(temp);
figure
imshow(temp)
figure
imshow(bw)
se = strel('disk',15);

mask_SHG = imdilate(bw,se);
figure, imshow(mask_SHG);

%%%%%%%%%%%%%%mask creations
mask_HE=im2bw(temp);

m=~mask_HE.*~mask_SHG;
mask2=m+mask_SHG;

mask2 = ~(bwareaopen(~mask2, 1000));

rgbImage=a_HE;
mask2 = bwconvhull(~mask2,'objects');
mask2 = bwconvhull(mask2,'objects');
% bw1=mask_HE&mask_SHG;
figure, imshow(mask2);
%
% imshow(im2bw(temp).*bw1)

mask=uint8(~mask2);
redPlane = rgbImage(:, :, 1);
greenPlane = rgbImage(:, :, 2);
bluePlane = rgbImage(:, :, 3);

% Do the masking.
maskedRed = redPlane .* mask;
maskedGreen = greenPlane .* mask;
maskedBlue = bluePlane .* mask;

% Combine back into a masked RGB image.
maskedRgbImage = cat(3, maskedRed, maskedGreen, maskedBlue);

figure, imshow(maskedRgbImage)

bw1=bwareaopen(bw,200);
figure
imshow(bw1+mask2)




st = regionprops(mask2, 'BoundingBox' );
mkdir('SHG subimages')
figure;
imshow(mask2)
for k = 1 : length(st)
    thisBB = st(k).BoundingBox;
    rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
        'EdgeColor','r','LineWidth',2 );
end


%%%%%%%%%%%%%%%creating rectangle around box

[m n]=size(bw1);

for k = 1:length(st)
    
    rect_cor(k,:)=st(k).BoundingBox;
    
    up_left_x=round(rect_cor(k,1)-rect_cor(k,3)/2);
    if up_left_x<1
        up_left_x=1;
    end
    
    up_left_y=round(rect_cor(k,2)-rect_cor(k,4)/2);
    if up_left_y<1
        up_left_y=1;
    end
    
    low_right_x=round(rect_cor(k,1)+rect_cor(k,3)*1.5);
    
    if low_right_x>n
        low_right_x=n;
    end
    
    low_right_y=round(rect_cor(k,2)+rect_cor(k,4)*1.5);
    
    if low_right_y>m
        low_right_y=m;
    end
    
    global_upleft_x(k)=up_left_x;
    global_upleft_y(k)=up_left_y;
    
%     path=['C:\Users\Adib\Google Drive\pattern recognition project\For AdibSagar\SHG subimages\Image' num2str(k) '.tiff'];
    %field=['Image' num2str(k)]
    value=a_SHG(up_left_y:low_right_y,up_left_x:low_right_x);
    obj_img.(sprintf('Image_%d', k))=value;
    
    I = mat2gray(obj_img.(sprintf('Image_%d', k)));
%     imwrite(I,path);
    
    figure
    imshow(I)
    
    [H_subimage(k).subimage, theta(k).subimage, rho(k).subimage] = hough(I);
end


%%block generaiotn code
% start_coor_x= global_upleft_x;
% start_coor_y= global_upleft_y;

for i=1:k
    
    [row(i) col(i)]=size(obj_img.(sprintf('Image_%d', i)));
    
    block_step_x=floor(col(i)/4);
    block_step_y=floor(row(i)/4);
    
    
    for j=1:5
        block_coord_x(i,j)=(j-1)*block_step_x+global_upleft_x(i);
        block_coord_y(i,j)=(j-1)*block_step_y+global_upleft_y(i);
        if block_coord_x(i,j)>n
            block_coord_x(i,j)=n;
        end
        if block_coord_y(i,j)>m
            block_coord_y(i,j)=m;
        end
    end
    
    
    
    
    
end


count=1;


for i=1:k
    x_coord=block_coord_x(i,:);
    y_coord=block_coord_y(i,:);
    for j=1:3
        for l=1:3
            block_coords(count,:)=[x_coord(j) x_coord(j+2) y_coord(l) y_coord(l+2)];
            
            block_association(count,:)=[i  x_coord(j) x_coord(j+2) y_coord(l) y_coord(l+2) j l];
            count=count+1;
        end
        
        
    end
end


%%%hough transform on SHG_images

for i=1:length(block_association)
    
    
    block(i).image=a_SHG(block_association(i,4):block_association(i,5),block_association(i,2):block_association(i,3));
    [H(i).hough, theta(i).hough, rho(i).hough] = hough(block(i).image);
    
end

for i=1:length(block_association)
    
    %     hough_max(i)=max(max((H(i).hough)));
    hough_feature.hough_max_SHG(i,1)=max(max((H(i).hough)));
    hough_feature.hough_avg_SHG(i,1)=mean(mean((H(i).hough)));
    hough_feature.hough_nonzero_avg_SHG(i,1)=mean(H(i).hough(find(H(i).hough~=0)));
    
end

A=[hough_feature.hough_max_SHG hough_feature.hough_avg_SHG hough_feature.hough_nonzero_avg_SHG];




 end






