% SVD example
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

[u,d,v]=svd(A);
n = size(A, 2);

scatter3(A(1,:),A(2,:),A(3,:),'.')

hold on
pc1 = d(1,1)/sqrt(n)*[0 0 0;u(1,1) u(2,1) u(3,1)];
t=plot3(pc1(:,1),pc1(:,2),pc1(:,3),'m');
set(t,'linewidth',4)
pc2 = d(2,2)/sqrt(n)*[0 0 0;u(1,2) u(2,2) u(3,2)];
t=plot3(pc2(:,1),pc2(:,2),pc2(:,3),'g');
set(t,'linewidth',4)
alpha(0.1)


