function [ bhat ] = ls( A, b )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

bhat = A * inv(transpose(A) * A) * transpose(A) * b;

end

