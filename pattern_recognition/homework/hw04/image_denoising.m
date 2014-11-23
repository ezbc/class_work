
clear
close all
n = 64;

% noise free image
x = double(phantom(n))*256;

% add noise
y = x+randn(size(x))*15;

% denoise by non-local means
sig = 150; % weighing factor; bigger sig = more smoothing
h = 5; % patch sidelength
x_nlm = nlm(y,h,sig);


% denoise by distance-based averaging
w = [1/16 1/8 1/16;1/8 1/4 1/8;1/16 1/8 1/16];
for i=1:n
    for j=1:n
        if (i==1)|(i==n)|(j==1)|(j==n)
            xavg(i,j) = y(i,j);  % don't process pixels at edge
        else
            b = [y(i-1,j-1) y(i,j-1) y(i+1,j-1); 
                 y(i-1,j) y(i,j) y(i+1,j);
                 y(i-1,j+1) y(i,j+1) y(i+1,j+1)];
            xavg(i,j) = sum(sum(b.*w));
        end
    end
end


figure(1);clf;
subplot(221);imagesc(x,[0,256]);axis image;colormap gray; title('original')
subplot(222);imagesc(y,[0,256]);axis image;colormap gray; title('noisy')
subplot(223);imagesc(xavg,[0,256]);axis image;colormap gray; title('averaging')
subplot(224);imagesc(x_nlm,[0,256]);axis image;colormap gray; title('non-local means')
linkaxes
