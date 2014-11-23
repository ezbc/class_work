A1 = [3 1 2; 0 3 3; 0 4 4; 6 1 4];
A2 = [1 1 2; 0 3 3; 0 4 4; 3 1 4];
A3 = randn(10, 6);

gs_ortho(A1)
rank(A1)
gs_ortho(A2)
rank(A2)
gs_ortho(A3)
rank(A3)


