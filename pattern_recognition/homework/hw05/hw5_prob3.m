clear
close all

data_ratings = load('jesterdata.mat');
data_new = load('newuser.mat');

% Extract data from structs
A = data_ratings.X;
b = data_new.b;
b_true = data_new.trueb;

resids = zeros(size(A, 2), 1);

for i = 1:size(A, 2);
    resids(i) = sum((A(b>-99, i) - b(b>-99)).^2);
end

'Index minimum residual'
[resids_sorted, indices] = sort(resids,'descend');

indices(1)
indices(2)



