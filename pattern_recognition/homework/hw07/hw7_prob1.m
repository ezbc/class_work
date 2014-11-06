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

% perform duel least squares
lambda = 10^-5;
xLS = A' * inv(A*A' + lambda*eye(size(A,1))) * b;
bLS = sign(A*xLS);

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

% Plot least squares
subplot(122); hold on;
for i=1:m
    a = A(i,:);
    if bLS(i)==1
        plot(a(1),a(2),'b.');
    else
        plot(a(1),a(2),'r.');
    end
end
axis('square')
title('least squares')



