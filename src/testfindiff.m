% Testing script for the findiffapprox function
% Implementation - C. W. Royer, January 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

A = [1000 -1 0;-1 2 -1;0 -1 2];
myfun = @(x) 0.5*x'*A*x;
%
for h=1:10
    fprintf('h=10^{-%d}\n',10^(-h));
    [H,flag] = findiffapprox(zeros(3,1),myfun,h);
    fprintf('Error=%1.2e\n',norm(H-A));
end
