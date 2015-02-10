clear;
close all;

load('wiscnet.mat')

E = mean(X, 2);

C = cov(X');

[V, D] = eig(C);

variability = 0.0;

% Total variability is the sum of the eigenvalues
tot_variability = sum(sum(D));
i = length(D) - 1;

while variability < 0.95 * tot_variability

    variability = sum(sum(D(i:length(D), i:length(D))));

    i = i - 1;

end



