clear
close all

data_ratings = load('jesterdata.mat');
data_new = load('newuser.mat');

% Extract data from structs
A = data_ratings.X;
b = data_new.b;
b_true = data_new.trueb;
A_sub = A(:, 1:100);

% Compute p with the power method
p = power_iter(A_sub);

% Compute the SVD, excluding rows with zeros, i.e, economy
[U, S, V] = svd(A_sub, 0);

% 
'1-norm between 1st col of U and p'
norm(U(:, 1) - p, 1)
'1-norm between 1st col of V and p'
norm(V(:, 1) - p, 1)




