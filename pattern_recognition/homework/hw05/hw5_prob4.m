clear
close all

data_ratings = load('jesterdata.mat');
data_new = load('newuser.mat');

% Extract data from structs
A = data_ratings.X;
b = data_new.b;
b_true = data_new.trueb;

% First attempt with SVD failed.
% ------------------------------
% Compute the SVD, excluding rows with zeros, i.e, economy
[U, S, V] = svd(A, 0);

% Plot the spectrum of X
fig = figure(1);clf;
subplot(111); scatter(linspace(1,size(S, 1), size(S, 1)), diag(S)); hold on;
xlabel('sigma')
ylabel('Value')
hold off;
%legend('data', 'eps = 0.001', 'eps = 0.01', 'eps = 0.1', 'Location', 'North')
box on;

saveas(fig, 'prob4_fig', 'png')



