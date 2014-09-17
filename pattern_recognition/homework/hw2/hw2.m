clear
close all
n = 64;

% noise free image
original= double(phantom(n)) * 256;

% add noise
noisey = original+ randn(size(original)) * 15;

% denoise by distance-based averaging
w = [1/16 1/8 1/16;1/8 1/4 1/8; 1/16 1/8 1/16];

for i=1:n
    for j=1:n
        if (i==1)||(i==n)||(j==1)||(j==n)
            xavg(i,j) = noisey(i,j);
        else
            b = [noisey(i-1, j-1) noisey(i, j-1) noisey(i+1, j-1);
                 noisey(i-1, j) noisey(i, j) noisey(i+1, j);
                 noisey(i-1, j+1) noisey(i, j+1) noisey(i+1, j+1)];
            xavg(i, j) = sum(sum(b.*w));
        end
    end
end

% Denoise by intensity similarity
intensity_list = reshape(noisey, numel(noisey), 1);
[intensity_list, sort_indices] = sort(intensity_list);
denoised_image = zeros(numel(intensity_list), 1);

for i=1:numel(intensity_list)
    int_val = intensity_list(i);
    
    diff = intensity_list - int_val;
    [diff, diff_indices] = sort(abs(diff));
    int_val_near = mean(intensity_list(diff_indices(2:1000)));

    int_val_avg = (int_val + int_val_near) / 2.0;

    denoised_image(sort_indices(i)) = int_val_avg;
end

shape = size(noisey);
denoised_image = reshape(denoised_image, shape(1), shape(2));

figure(1);clf;
subplot(141);imagesc(original, [0,256]);axis image; colormap gray
subplot(142);imagesc(noisey, [0,256]);axis image; colormap gray
subplot(143);imagesc(xavg, [0,256]);axis image; colormap gray
subplot(144);imagesc(denoised_image, [0,256]);axis image; colormap gray
linkaxes


