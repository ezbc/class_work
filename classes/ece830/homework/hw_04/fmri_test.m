

% THIS FILE CAN BE DOWNLOADED AT http://www.ece.wisc.edu/~nowak/
% fMRI detection
% R. Nowak
% 10/23/08

clear
close all

% load fMRI data and reference signal waveform

% loads 122 x 64 x 64 array called 'x', consisting of 122 brain images
% each image is 64x64 pixels 
load fmri.mat

% loads 122 x 1 reference signal waveform 
load ref.mat

figure(1)
plot(s)
axis([1,122,0,1500])
title('reference signal: time series of very active voxel')

% remove mean (nuissance parameter) from reference signal 
s = s - mean(s);

% compute "correlation" image
% correlate each pixel time-series with reference signal
% this image is an array of test statistics
%figure(2)
%for i=1:122
%    im(1:64,1:64) = x(i,:,:);
%    imagesc(im);
%    axis('square')
%    colormap('gray')
%    title(['Time ',num2str(i*2),' seconds'])
%    pause(.25);
%end
%
t = zeros(64,64);
for i=1:64
for j=1:64
	t(i,j) = sum(x(:,i,j).*s);
    mx(i,j) = sum(x(:,i,j));
end
end
mx = mx/length(s);

% display correlation image
figure(3)
imagesc(t)
axis('square')
colormap('gray')
title('correlations of all voxels with reference')

% compare test statistics to threshold (set at 35% of max value)
% array d consists of 0's and 1's, corresponding to 
% decide H_0 = no activation, decide H_1 = activation
gam = 0.35*max(max(t));
d = (t>gam);

% display detection results
% pseudo-color red indicates activity
% pseudo-color blue (weighted by mean pixel values) indicates regions
% with no activation
figure(4)
 
dmap(:,:,1) = d .* (mx>400);
dmap(:,:,2) = (1-d).*(mx>400)/2;
dmap(:,:,3) = (1-d).*(mx>400)/2;
imagesc(dmap)
axis('square')
title('crude map of activations')


% Look at uncorrelated voxels
t = zeros(64,64);
for i=1:64
for j=1:64
	t(i,j) = sum(x(:,i,j).*s);
    mx(i,j) = sum(x(:,i,j));
end
end
mx = mx/length(s);

% display correlation image
figure(4)
plot(x(:,10,10))
axis('square')
colormap('gray')
title('noisy voxel')

