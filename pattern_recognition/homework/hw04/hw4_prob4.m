data = load('face_emotion_data.mat');

A = data.A;
b = data.b;

Q = gs_ortho(A);
bhat_orth = Q * transpose(Q) * b;
bhat_ls = ls(A, b);
'1a'
'The Pearson Correlation coefficient between the two solutions is '
p = corr(bhat_orth, bhat_ls)

'1b'
O = orth(A); 
bhat_orth = O * transpose(O) * b;
bhat_gsortho = Q * transpose(Q) * b;

'The Pearson Correlation coefficient between the two solutions is '
p = corr(bhat_gsortho, bhat_orth)

