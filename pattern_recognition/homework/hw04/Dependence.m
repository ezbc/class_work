%This program determines whether a set of column vectors is linearly 
%independent
%or linearly dependent. It accepts a Matrix
%"B" and returns a scalar "d" which equals "1" if
%the columns of "A" are Linearly Independent and "0" if they are 
%Linearly Dependent.
function [d]=Dependence(B)
C=rref(B);
m=length(diag(B(:,1)));
n=length(B(1,:));
if n>m
 d=0;
else
s=sum(diag(C));
if n>s
d=0;
else
d=1;
end
end