clear
close all

data_ratings = load('jesterdata.mat');
data_new = load('newuser.mat');

% Extract data from structs
A = data_ratings.X;
b = data_new.b;
b_true = data_new.trueb;
A_sub = A(:, 1:100);

b_i = A(:);

% Compute p with the power method
p = power_iter(A_sub, b_i);

% Compute the SVD, excluding rows with zeros, i.e, economy
[U, S, V] = svd(A_sub, 0);






