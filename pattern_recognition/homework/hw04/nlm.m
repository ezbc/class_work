
% nonlocal means
% rebecca willett, 9/15/2014

function y = nlm(x,h,sig)

h = floor(h/2); % h must be odd

[n,m] = size(x);
x_pad = padarray(x,[h h],'replicate','both');

y = zeros(size(x));

sig = sig*(2*h+1)^2;

fprintf(' i = ')
for i = 1:n
  for j = 1:m
    patch_ij = x_pad((i):(i+2*h),(j):(j+2*h));
    patch_sum = zeros(2*h+1);
    weight_sum = 0;
    for k = 1:n
      for l = 1:n
        patch_kl = x_pad((k):(k+2*h),(l):(l+2*h));
        dist = norm(patch_ij-patch_kl,'fro')^2;
        weight_kl = exp(-dist/sig);
        weight_sum = weight_sum + weight_kl;
        patch_sum = patch_sum + patch_kl*weight_kl;
      end
    end
    y(i,j) = patch_sum(h+1,h+1)/weight_sum;
  end
  fprintf('%d... ',i);
  if mod(i,16)==0
    fprintf('\n i = ');
  end
end
fprintf('\n');


