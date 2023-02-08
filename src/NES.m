function [EstString,negativefound] = NES(H,orderoption,buildstyle,verboseflag)
%% NES (Negative Eigenvalue Search) by Warren Hare, 2021
% Input: (H, orderoption,buildstyle,verboseflag)
%         H = Matrix that may or may not have an negative eigenvalue
%         orderoption = flag that controls how the search proceeds
%         buildstyle = flag that controls how the build proceeds
%         verboseflag = sets visual output level to 
%                  very verbose(2), verbose (1), or quiet (0)
%
% Overview: The code begins by making the diagonal of H, D
%         It then expands entries in D->Hestimate, slowing filling it out 
%         to H.  The code stops when it feels that D+ provides evidence 
%         that H has negative eigenvalues, or provides evidence that H has 
%         no negative eigevalues.  Worst case, (n^2-n)/2 iterations will
%         cause Hestimate = H, so full knowledge is obtained.
%         
% Output: (EstString, stopiteration)
%         EstString = vector of approximated eigenvalues
%         stopiteration = iteration when the program stopped (i.e.,
%         detected a negative eigenvalue or guessed that all eigenvalues
%         are positive)

%format long
%% Set up
    D=diag(H);n=length(D);Hestimate=diag(D);
%% Set order of growth for H
    if n==1 % for 1 x 1 matrices, orderoption must = 1.
        ordopt=1;else ordopt=orderoption;
    end
    Order=SetOrder(D,n,ordopt,buildstyle);
    
%% Grows the matrix and update the eigenvalue estimations
    Hestimate=diag(D);
    EstString=eig(Hestimate);
    counter=0; negativefound=0;
    if verboseflag>1
        fprintf('At iteration %d, the minimal estimated eigenvalue is %2.2f\n', counter, min(EstString));
    end
%   Added by Clement - handle the easy case of finding negative curvature
%   right away
    if min(EstString)>=0
%      
        for counter=1:size(Order,2)
            i=Order(1,counter);
            j=Order(2,counter);
            Hestimate(i,j)=H(i,j);
            Hestimate(j,i)=H(j,i);
%           Correction - Use the largest principle submatrix or submatrices
%           contained in Hestimate
            [Hps,nps,sps] = findsubmatrix(Hestimate,Order,counter);
            emin = Inf;
            for ips=1:nps
                EstString = eig(Hps{ips});
                emin = min(emin,min(EstString));
            end
            if verboseflag>1
                fprintf('At iteration %d, the minimal estimated eigenvalue is %2.2f\n', counter, emin);
                fprintf('This estimate was computed based on %d by %d submatrices\n',sps,sps);
            end
            if negativefound==0 && emin<0
                negativefound=counter;
            end
        end
    end
    
    if verboseflag
        fprintf('The first negative eigenvalue was detected at iteration %d\n', negativefound);
    end
end
