x1 = [1; 0; 0;];
x2 = [1; 0.6; 0.8;];
Q = [x1 x2];
b = [1; 3; 1;];

bhat = Q * inv(Q' * Q) * Q' * b

