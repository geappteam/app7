function [ E ] = evalErrorPower( reference, prediction )
%EVALERRORPOWER Summary of this function goes here
%   Detailed explanation goes here

e = reference - prediction;

E = sum(e.^2);

end

