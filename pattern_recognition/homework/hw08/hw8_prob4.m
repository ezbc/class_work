clear;
close all;

% mean error of three methods
% LS 48.3%
% Poly kernel 3.93%
% Gauss kernel 1.35%

% generate data
m_list = [10, 100, 1000];
m_list = [10, 100, 1000];
b_list = {};
A_list = {};
n = 2;

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
        
        % check to make sure not all the same
        if b == b(i)
            b(1) = 0;
        end

        % LS
        % --
        xLS = A' * inv(A*A' + lambda*eye(size(A,1))) * b;
        b_hat_ls = sign(A*xLS);
        err = sum(b_hat_ls ~= b)/m;
        err_array(j, 1) = err;

        % Gauss kernel
        % ------------
        svm = svmtrain(A, b, 'kernel_function','rbf');
        b_hat_gauss = svmclassify(svm, A);

        err = sum(b_hat_gauss ~= b)/m;
        err_array(j, 2) = err;

        % Poly kernel
        % ------------
        svm = svmtrain(A, b, 'kernel_function','polynomial');

        xhat = svm.Alpha;
        b_hat_poly = svmclassify(svm, A);

        err = sum(b_hat_poly ~= b)/m;
        err_array(j, 3) = err;
    end
    err_list{i} = err_array;
end

save('err_list', 'err_list');

text = {'LS', 'Poly', 'Gauss'};
for i=1:size(m_list,2)
    for j=1:3
        disp(m_list(i))
        disp(text{j})
        s = sum(err_list{i});
        disp(s(j))
    end
end

% Plot data
figure(1); 
subplot(141); hold on;
for i=1:m
    a = A(i, :);
    if b(i)==1
        plot(a(1),a(2),'b.');
    else
        plot(a(1),a(2),'r.');
    end
end
axis('square')
title('Training Data')

subplot(142); hold on;
for i=1:m
    a = A(i,:);
    if b_hat_ls(i)==1
        plot(a(1),a(2),'b.');
    else
        plot(a(1),a(2),'r.');
    end
end
axis('square')
title('Least Squares')

% Plot solution
subplot(143); hold on;
for i=1:m
    a = A(i,:);
    if b_hat_poly(i)==1
        plot(a(1),a(2),'b.');
    else
        plot(a(1),a(2),'r.');
    end
end
axis('square')
title('Polynomial Kernel')

subplot(144); hold on;
for i=1:m
    a = A(i,:);
    if b_hat_gauss(i)==1
        plot(a(1),a(2),'b.');
    else
        plot(a(1),a(2),'r.');
    end
end
axis('square')
title('Gaussian Kernel')


