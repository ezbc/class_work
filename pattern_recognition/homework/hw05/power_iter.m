% Calculates eigenvalue of matrix A using power iteration method
function b = power_iter(A, b_i)

    b_i = zeros(size(A, 1), 1);
    b_k = ones(size(A, 1), 1);

    toler = 1e-10;

    while norm(b_i - b_k, 1) > toler
        b_i = b_k;
        b_k = A * b_i / norm(A * b_i);
    end

    b = b_k;
end
