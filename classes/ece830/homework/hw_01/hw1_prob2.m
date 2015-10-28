clear;
close all;

n = 32; B = 1/16; M = 1000;
S = zeros(n,n);
for m=1:M;
    X = cos(2*pi*(2*rand-1)*B*(0:n-1)'+pi*rand);
    Y = X + 0.5 * randn(n,1);
    S = S+Y*Y';
end;

S = S/M;

lambda = eig(S); 

disp('Number of dimensions in eigenspace:')
disp(length(lambda(find(lambda > 0.5))))



