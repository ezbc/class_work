clear
close all

data_ratings = load('jesterdata.mat');
data_new = load('newuser.mat');

% Extract data from structs
A = data_ratings.X;
b = data_new.b;
b_true = data_new.trueb;
A_sub = A(:, 1:20);

% First attempt with SVD failed.
% ------------------------------
% Compute the SVD, excluding rows with zeros, i.e, economy
[U, S, V] = svd(A_sub, 0);

% Use least squares solution for overdetermined solution
x_hat = U * (S\(V(b > -99, :)' * b(b > -99)));

% Get the best estimate of b
disp('SVD best estimate')
b_hat = A_sub * x_hat;

% Measure differences between true and measured b
diff = b_true - b_hat;
disp('Difference in maximum location of b_true')
disp(diff(29))

disp('std difference')
disp(std(diff))

disp('Maximum location of b_hat')
disp(find(b_hat==max(b_hat)))

% Plot the figure
fig = figure(1);clf;
subplot(111); scatter(linspace(1,100),diff); hold on;

xlabel('m')
ylabel('diff')
hold off;
%legend('data', 'eps = 0.001', 'eps = 0.01', 'eps = 0.1', 'Location', 'North')
box on;

saveas(fig, 'prob1_fig', 'png')





