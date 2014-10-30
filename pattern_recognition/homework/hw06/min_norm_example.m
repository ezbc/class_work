% min norm example
clear
close all

load face_emotion_data

% generate k new features (random combo of other features)
k=3;
a = A*randn(9,k);

% new data matrix
B = [A a];

% svds
svd(A)
svd(B)


% LS solutions
b1 = A*inv(A'*A)*A'*b;

b2 = B*inv(B'*B)*B'*b;

% min norm solution
[U,D,V] = svd(B);
Dinv = zeros(size(D));
for i=1:min(size(D))
    if D(i,i)>10e-10,
        Dinv(i,i)=1/D(i,i);
    end
end
b3 = B*V*Dinv'*U'*b;

figure(1)
subplot(121)
imagesc(D)
axis equal
subplot(122)
imagesc(Dinv')
axis equal

b3 = B*pinv(B)*b;

% regularized solution
lam = 1;
b4 = B*inv(B'*B+lam*eye(size(B'*B)))*B'*b;


figure(2)
subplot(411)
plot(b)
subplot(412)
plot(b1)
ax= axis;
subplot(413)
plot(b2)
axis(ax)
subplot(414)
plot(b3)
axis(ax)

errs = [sum(b~=sign(b1)) sum(b~=sign(b2)) sum(b~=sign(b3)) sum(b~=sign(b4))]
