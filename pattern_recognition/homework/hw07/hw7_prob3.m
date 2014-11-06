clear;
close all;

% generate data
m = 1000;
n = 2;
b = zeros(m,1);
for i=1:m
    a = 2*rand(2,1)-1;
    A(i,:)=a';
    b(i) = sign(a(1)^2+a(2)^2-.5);
end

% perform polynomial kernel estimation of b
lambda = 10^-5;
K = zeros(m, m);
for i = 1:m
    for j = 1:m
        K(i, j) = (dot(A(i, :)', A(j, :)) + 1)^2;
    end
end

alpha = (K + lambda * eye(m)) \ b;
x_hat = A' * alpha;
b_hat = sign(K*alpha);

% Plot data
figure(1); 
subplot(121); hold on;
for i=1:m
    a = A(i, :);
    if b(i)==1
        plot(a(1),a(2),'b.');
    else
        plot(a(1),a(2),'r.');
    end
end
axis('square')
title('training data')

% Plot solution
subplot(122); hold on;
for i=1:m
    a = A(i,:);
    if b_hat(i)==1
        plot(a(1),a(2),'b.');
    else
        plot(a(1),a(2),'r.');
    end
end
axis('square')
title('Polynomial Kernel')



