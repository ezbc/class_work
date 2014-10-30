clear
close all

data = load('face_emotion_data.mat');
A = data.A;
b = data.b;

% Number of partitions
npar = 8;
N = size(A, 1);

% Scatter row positions
pos   = randperm(N);

% Bin the positions into npar partitions
edges = round(linspace(1, N+1, npar+1));

%Now you can "physically" partition A, or apply your code to the segments of
%without actually separating into blocks.

% Partition A
prtA  = cell(npar,1);
prtb  = cell(npar,1);
for ii = 1:npar
    idx      = edges(ii):edges(ii+1)-1;
    prtA{ii} = A(pos(idx), :); % or apply code to the selection of A
    prtb{ii} = b(pos(idx)); % or apply code to the selection of A
end

% Grab unique set indices for partitions
indices = zeros(56, 2);
count = 1;
for test_index = 1:8
    for tune_index = 1:8
        if test_index ~= tune_index
            indices(count, :) = [test_index tune_index];
            count = count + 1;
        end
    end
end

% Define list of regularization parameter lambda
lambdas = [0, 2^-1 2^0 2^1 2^2 2^3 2^4];
K = length(lambdas);

% For each partition, estimate x_hat with SVD
errs = zeros(size(indices, 1), 1);
lambda_mins = zeros(size(indices, 1), 1);
for i = 1:size(indices, 1)
    % Snag test and tuning dataset
    A_test = prtA{indices(i, 1)};
    A_tune = prtA{indices(i, 2)};
    b_test = prtb{indices(i, 1)};
    b_tune = prtb{indices(i, 2)};

    % for indices not in the test + tune indices, create train
    train_indices = setdiff(linspace(1, 8, 8), indices(i, :));
    A_train = cat(1, prtA{train_indices});
    b_train = cat(1, prtb{train_indices});

    % Create identity matrix for subset
    I_train = eye(size(A_train, 2));

    % Calculate least squares solution for each regularization
    x_hats = cell(K, 1);
    for k = 1:K
        lambda = lambdas(k);
        x_hats{k} = inv(A_train.' * A_train + lambda * I_train) * ...
                    A_train.' * b_train;
    end

    % Tune the regularization parameter
    resid_min = Inf;
    lambda_min = 1;
    min_index = 1;
    for k = 1:K
        resid = norm(b_tune - A_tune * x_hats{k});
        if resid < resid_min
            resid_min = resid;
            lambda_min = lambdas(k);
            min_index = k;
        end
    end
    lambda_mins(i) = lambda_min;

    % Test tuned model on test data
    errs(i) = norm(b_test - A_test * x_hats{min_index});
end

% calculate mean of errors
%disp('Lambda_mins')
%disp(lambda_mins);
disp('Mean of errors')
disp(mean(errs))








