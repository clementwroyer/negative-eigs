function [HFD,flag]=findiffapprox(x,fun,h)
%
% A basic code to compute a finite-difference estimates of a Hessian matrix.
% Uses a centered finite-difference formula based on coordinate vectors.
%
% Inputs:
%  x: Point at which the Hessian is to be estimated
%  fun: Function handle of a twice differentiable function
%  h: Stepsize for the finite difference formulas.
%
% Output:
%  HFD: Finite-difference matrix
%  flag: Success of the method
%   0 - Estimate was successfully computed
%   m>0 - m infinite or NaN value were obtained
%
% Implementation: C. W. Royer, January 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(x);
I = eye(n);
HFD = zeros(n);
flag = 0;
%
fx = feval(fun,x);
fpm = zeros(n,2);
% Diagonal part
for i=1:n
    fpm(i,1) = feval(fun,x+h*I(:,i));
    if isinf(fpm(i,1)) || isnan(fpm(i,1))
        flag = flag+1;
    end
    fpm(i,2) = feval(fun,x+2*h*I(:,i));
%    fpm(i,2) = feval(fun,x-h*I(:,i));
    if isinf(fpm(i,2)) || isnan(fpm(i,2))
        flag = flag+1;
    end
%    HFD(i,i) = (fpm(i,1)-2*fx+fpm(i,2))/(h^2);
    HFD(i,i) = (fpm(i,2)-2*fpm(i,1)+fx)/(h^2);
end
% Off-diagonal terms
for i=1:n
    for j=(i+1):n
        ft = feval(fun,x+h*I(:,i)+h*I(:,j));
        if isinf(ft) || isnan(ft)
            flag = flag+1;
        end
        HFD(i,j) = (ft-fpm(i,1)-fpm(j,1)+fx)/(h^2);
        HFD(j,i) = HFD(i,j);
    end
end     
% Symmetrization
HFD = 0.5*(HFD+HFD');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
