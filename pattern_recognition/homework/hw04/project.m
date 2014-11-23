function [ w ] = project( u, v)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

w = (tranpose(u) * v) / (transpose(u) * u) * u;

end

