clear;
close all;

% mean error of three methods
% LS 48.3%
% Poly kernel 3.93%
% Gauss kernel 1.35%

% generate data
m_list = [10, 100, 1000];
b_list = {};
A_list = {};
n = 2;
for j = 1:size(m_list,2)
    m = m_list(j);
    b = zeros(m,1);
    for i=1:m
        a = 2*rand(2,1)-1;
        A(i,:)=a';
        b(i) = sign(a(1)^2+a(2)^2-.5);
    end
    b_list{j} = b;
    A_list{j} = A;
end

% perform polynomial kernel estimation of b
lambda = 10^-5;
n_runs = 100;
err_list = {};
for i = 1:size(m_list, 2)
    m = m_list(i);
    err_array = zeros(n_runs,3);
    for j = 1:n_runs

        % Generate data
        b = zeros(m,1);
        A = zeros(m, 2);
        for ii=1:m
            a = 2*rand(2,1)-1;
            A(ii,:) = a';
            b(ii) = sign(a(1)^2+a(2)^2-.5);
        end

        % LS
        % --
        xLS = A' * inv(A*A' + lambda*eye(size(A,1))) * b;
        b_hat = sign(A*xLS);
        err = sum(b_hat ~= b)/m;
        err_array(j, 1) = err;

        % Polynomial kernel
        % -----------------
        K = zeros(m, m);
        for ii = 1:m
            for jj = 1:m
                K(ii, jj) = (dot(A(ii, :)', A(jj, :)) + 1)^2;
            end
        end
        alpha = (K + lambda * eye(m)) \ b;
        x_hat = A' * alpha;
        b_hat = sign(K*alpha);
        err = sum(b_hat ~= b)/m;
        err_array(j, 2) = err;

        % Gaussian kernel
        % ---------------
        K = zeros(m, m);
        for ii = 1:m
            for jj = 1:m
                K(ii, jj) = exp(-0.5 * norm(A(ii, :) - A(jj, :))^2);
            end
        end

        alpha = (K + lambda * eye(m)) \ b;
        x_hat = A' * alpha;
        b_hat = sign(K*alpha);
        err = sum(b_hat ~= b)/m;
        err_array(j, 3) = err;
    end
    err_list{i} = err_array;
end

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


