%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing script - Detecting negative curvature in Hessian matrices from 
% CUTEst.
%
% This script operates on matrices from CUTEst that possess negative 
% eigenvalues. For each matrix, the goal is to detect negative curvature using 
% submatrices. All possible orderings for selecting the coefficients defining 
% the submatrices are compared.
%
% Started August 31, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
close all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
findiff=1;
% Load Hessian matrices
if findiff
    load('HESSIANSNCFD')
else
    load('HESSIANSNC')
end
npbs = length(pbdims);
nmats = npbs;
% Set dimensions 
%maxdim = 50;
maxdim = 8;
mindim = 2;
% Select problems with appropriate dimensions
Imat = find(pbdims<=maxdim & pbdims>=mindim);
% Random
randorthog=0;
%randorthog=1;
% Select the ordering options (avoid combinatorial explosion)
if maxdim>8
    nocombi=1;
else
    nocombi=0;
end
%nocombi=1;
verbose=2;%Verbose level in FindBestOrder subroutine
% Plot all combinations in combinatorial case?
plotcombi=0;
nsel = length(pbdims(Imat)); % Number of problems selected
if findiff 
    nFD = length(hFD);
else
    nFD = 0;
end
% Total number of matrices (original+Newton iterations)
nmat = nsel*(1+nitsN)*(1+nFD);
% TO DO - Adjust numbering
%
% Pre-allocate output data structures
%   minOrd* represents the minimal value of a certain ordering
%   bestOrd* contains the actual ordering that achieves the best value
minOrd1 = zeros(nmat,1);
bestOrd1 = cell(nmat,1);
minOrd2 = zeros(nmat,1);
bestOrd2 = cell(nmat,1);
%
Ikeep = [];% Problem indices
%%%%%%%%%%%%%%%%%%%
% Main loop
for i=1:nsel
    auxi = (nFD+1)*(nitsN+1)*(i-1)+1;
    fprintf('Problem %s\n',pbnames{Imat(i)});
    fprintf('\t Init pt (Exact) - index %d \n',auxi);
    myH = pbmats{Imat(i)};
    if randorthog
        myv = randn(pbdims(Imat(i)),1);
        [myQ,~] = qr(myv);
        myH = myQ*myH*myQ';
    end 
%   Compute the best orderings
%    [minOrd1(auxi),bestOrd1{auxi},minOrd2(auxi),...
%    bestOrd2{auxi}] = FindBestOrder(pbmats{Imat(i)},verbose,nocombi);
    [minOrd1(auxi),bestOrd1{auxi},minOrd2(auxi),...
    bestOrd2{auxi}] = FindBestOrder(myH,verbose,nocombi);
    Ikeep = [Ikeep auxi];
    for k=1:nFD
        fprintf('\t Init Pt (FD=%1.0e) - index %d \n',hFD(k),auxi+k);
        if pbeigsFD(Imat(i),k)<0
            myH = pbmatsFD{Imat(i)}{k}
            if randorthog
                myv = randn(pbdims(Imat(i)),1);
                [myQ,~] = qr(myv);
                myH = myQ*myH*myQ';
            end 
            [minOrd1(auxi+k),bestOrd1{auxi+k},minOrd2(auxi+k),...
            bestOrd2{auxi+k}] = FindBestOrder(myH,verbose,nocombi);
            Ikeep = [Ikeep auxi+k];
        else
            minOrd1(auxi+k)=-1;
            bestOrd1{auxi+k}=[];
            minOrd2(auxi+k)=-1;
            bestOrd2{auxi+k}=[];
        end
    end
    auxi = auxi+nFD;
    for j=1:nitsN
        auxj = auxi+(nFD+1)*(j-1)+1;
        fprintf('\t It Newton %d - index %d \n',j,auxj);
        if pbeigsN(Imat(i),j)<0
            myH = pbmatsN{Imat(i)}{j};
            if randorthog
                myv = randn(pbdims(Imat(i)),1);
                [myQ,~] = qr(myv);
                myH = myQ*myH*myQ';
            end 
%            [minOrd1(auxi+j),bestOrd1{auxi+j},minOrd2(auxi+j),...
%            bestOrd2{auxi+j}] = FindBestOrder(...
%            pbmatsN{Imat(i)}{j},verbose,nocombi);
            [minOrd1(auxj),bestOrd1{auxj},minOrd2(auxj),...
            bestOrd2{auxj}] = FindBestOrder(myH,verbose,nocombi);
            Ikeep = [Ikeep auxj];
        else
            minOrd1(auxj)=-1;
            bestOrd1{auxj}=[];
            minOrd2(auxj)=-1;
            bestOrd2{auxj}=[];
        end
        for k=1:nFD
            fprintf('\t It Newton %d (FD=%1.0e) -index %d\n',j,hFD(k),auxj+k); 
            if pbeigsNFD(Imat(i),j,k)<0
                myH = pbmatsNFD{Imat(i)}{j}{k};
                if randorthog
                    myv = randn(pbdims(Imat(i)),1);
                    [myQ,~] = qr(myv);
                    myH = myQ*myH*myQ';
                end 
                [minOrd1(auxj+k),bestOrd1{auxj+k},minOrd2(auxj+k),...
                bestOrd2{auxj+k}] = FindBestOrder(myH,verbose,nocombi);
                Ikeep = [Ikeep auxj+k];
            else
                minOrd1(auxj+k)=-1;
                bestOrd1{auxj+k}=[];
                minOrd2(auxj+k)=-1;
                bestOrd2{auxj+k}=[];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%
%
% Save the relevant data into a .mat file
save DATANESCUTEST pbdims nitsN nFD Imat Ikeep minOrd1 bestOrd1 minOrd2 bestOrd2
%
%%%%%%%%%%%%%%%%%%%
% Write the desired outputs in a data file
fid = fopen('ResultsHessianCUTEst','w');
fprintf(fid,'Minimum allowed dimension: %d\n\n',mindim);
fprintf(fid,'Maximum allowed dimension: %d\n\n',maxdim);
if nocombi
    fprintf(fid,'Orderings:\n');
    fprintf(fid,'\t 1 - 1:n\n');
    fprintf(fid,'\t 2 - Smallest diag to largest diag\n');
    fprintf(fid,'\t 3 - Largest diag to smallest diag\n');
    fprintf(fid,'\t 4 - Flipflop smallest/largest diag\n\n');
else
    fprintf(fid,'Orderings: Combinatorial\n\n');
end
% Determine winners between Build 1 and Build 2
nkeep = length(Ikeep);
nkeepm2 = 0;
nkeepm3 = 0;
bestmin1 = (minOrd1<minOrd2);
sb1 = sum(bestmin1(Ikeep));
bestmin2 = (minOrd1>minOrd2);
sb2 = sum(bestmin2(Ikeep));
fprintf(fid,...
'Build type 1 gave the best order for %d out of %d matrices.\n',...
sb1,nkeep);
fprintf(fid,...
'Build type 2 gave the best order for %d out of %d matrices.\n',...
sb2,nkeep);
fprintf(fid,...
'Both builds gave the best order for %d out of %d matrices.\n\n',...
nkeep-sb1-sb2,nkeep);
% Counters in the case of dimension-independent strategies
if nocombi
    countbest1=[0 0 0 0];
    countbest1m2 = [0 0 0 0];
    countbest1m3 = [0 0 0 0];
    countbest2=[0 0 0 0];
    countbest2m2 = [0 0 0 0];
    countbest2m3 = [0 0 0 0];
end
%
fprintf(fid,'Problem \t| Dim | It Newton | FinDiff | MinEig |');
fprintf(fid,' BestBuild | MinUpdate | BestOrder\n');
fprintf(fid,'-----------------------------------------------------\n\n');
%%%%%%%%%
% Detailed results for all problems
% Outer loop on CUTEst problems, inner loop on the various matrices 
% with respect to that problem.
for i=1:nsel
    if ~nocombi 
        Pi = perms(1:pbdims(Imat(i)));
    end
    if nocombi
        nways=4;
    else
        nways=factorial(pbdims(Imat(i)));
    end
    auxi = (nFD+1)*(nitsN+1)*(i-1);
    auxj = auxi + 1;
    for j=0:1:nitsN
        for jFD=0:1:nFD
            fprintf('Problem %s ItN %d FD %d - index %d\n',...
            pbnames{Imat(i)},j,jFD,auxj+jFD);
            fprintf(fid,'%s\t%d \t %d \t',pbnames{Imat(i)},...
            pbdims(Imat(i)),j);
            if jFD==0
                fprintf(fid,'Exact \t');
            else
                fprintf(fid,'%1.0e \t',hFD(jFD));
            end
            if j==0 
                if jFD==0
                    fprintf(fid,'%1.3e \t\t',pbeigs(Imat(i)));
                else
                    fprintf(fid,'%1.3e \t\t',pbeigsFD(Imat(i),jFD));
                end
            else
                if jFD==0
                    fprintf(fid,'%1.3e \t\t',pbeigsN(Imat(i),j));
                else
                    fprintf(fid,'%1.3e \t\t',pbeigsNFD(Imat(i),j,jFD));
                end
            end 
            if ismember(auxj+jFD,Ikeep)
                if pbdims(Imat(i))>=3
                    nkeepm2 = nkeepm2+1;
                    if pbdims(Imat(i))>3
                        nkeepm3 = nkeepm3+1;
                    end
                end
                if bestmin1(auxj+jFD)
                    fprintf(fid,'1 \t %d \t',minOrd1(auxj+jFD));
                    for k=1:length(bestOrd1{auxj+jFD})
                        if nocombi || (~nocombi & plotcombi)
                            fprintf(fid,'%d, ',bestOrd1{auxj+jFD}(k));
                            if nocombi
                                countbest1(bestOrd1{auxj+jFD}(k)) = countbest1(...
                                bestOrd1{auxj+jFD}(k))+1;
                                if pbdims(Imat(i))>=3
                                    countbest1m2(...
                                    bestOrd1{auxj+jFD}(k)) = countbest1m2(...
                                    bestOrd1{auxj+jFD}(k))+1;
                                    if pbdims(Imat(i))>3
                                        countbest1m3(...
                                        bestOrd1{auxj+jFD}(k)) = countbest1m3(...
                                        bestOrd1{auxj+jFD}(k))+1;
                                    end
                                end
                            end
                        else
                            fprintf(fid,'[ ');
                            for l=1:pbdims(Imat(i))
                                fprintf(fid,'%d ',Pi(bestOrd1{auxj+jFD}(k),l));
                            end
                            fprintf(fid,'],');
                        end
                    end                
                    fprintf(fid,' \n');
                else
                    if bestmin2(auxj+jFD)
                        fprintf(fid,'2 \t %d \t',minOrd2(auxj+jFD));
                        for k=1:length(bestOrd2{auxj+jFD})
                            if nocombi || (~nocombi & plotcombi)
                                fprintf(fid,'%d, ',bestOrd2{auxj+jFD}(k));
                                if nocombi
                                    countbest2(...
                                    bestOrd2{auxj+jFD}(k)) = countbest2(...
                                    bestOrd2{auxj+jFD}(k))+1;
                                    if pbdims(Imat(i))>=3
                                        countbest2m2(...
                                        bestOrd2{auxj+jFD}(k)) = countbest2m2(...
                                        bestOrd2{auxj+jFD}(k))+1;
                                        if pbdims(Imat(i))>3
                                            countbest2m3(...
                                            bestOrd2{auxj+jFD}(k)) = countbest2m3(...
                                            bestOrd2{auxj+jFD}(k))+1;
                                        end
                                    end
                                end
                            end
                        end
                        fprintf(fid,' \n');
                    else
                        fprintf(fid,'1+2 \t %d \t',minOrd1(auxj+jFD));
                        for k=1:length(bestOrd1{auxj+jFD})
                            if nocombi || (~nocombi & plotcombi)
                                fprintf(fid,'%d, ',bestOrd1{auxj+jFD}(k));
                                if nocombi
                                    countbest1(...
                                    bestOrd1{auxj+jFD}(k)) = countbest1(...
                                    bestOrd1{auxj+jFD}(k))+1;
                                    if pbdims(Imat(i))>=3
                                        countbest1m2(...
                                        bestOrd1{auxj+jFD}(k)) = countbest1m2(...
                                        bestOrd1{auxj+jFD}(k))+1;
                                        if pbdims(Imat(i))>3
                                            countbest1m3(...
                                            bestOrd1{auxj+jFD}(k)) = countbest1m3(...
                                            bestOrd1{auxj+jFD}(k))+1;
                                        end
                                    end
                                end
                            end
                        end
                        fprintf(fid,'+');
                        for k=1:length(bestOrd2{auxj+jFD})
                            if nocombi || (~nocombi & plotcombi)
                                fprintf(fid,'%d, ',bestOrd2{auxj+jFD}(k));
                                if nocombi
                                    countbest2(...
                                    bestOrd2{auxj+jFD}(k)) = countbest2(...
                                    bestOrd2{auxj+jFD}(k))+1;
                                    if pbdims(Imat(i))>=3
                                        countbest2m2(...
                                        bestOrd2{auxj+jFD}(k)) = countbest2m2(...
                                        bestOrd2{auxj+jFD}(k))+1;
                                        if pbdims(Imat(i))>3
                                            countbest2m3(...
                                            bestOrd2{auxj+jFD}(k)) = countbest2m3(...
                                            bestOrd2{auxj+jFD}(k))+1;
                                        end
                                    end
                                end
                            end
                        end
                        fprintf(fid,' \n');
                    end
                end
            else %~ismember(auxi+j,Ikeep)
                fprintf(fid,'- \t - \t -\n');
            end
        end
        auxj = auxj+nFD+1;
    end
end
%
%%%%%%%%%%%%%%%%%%%
% Final statistics:
% For dimension-independent strategies, plot the percentage of pbms for 
% which each strategy was the best
fprintf(fid,'=============================================================\n');
if nocombi
    countbest1 = (100/nkeep)*countbest1;
    countbest2 = (100/nkeep)*countbest2;
    fprintf(fid,'Percentage of problems with best strategies [1 2 3 4] ');
    fprintf(fid,'(over %d matrices):\n',nkeep);
    fprintf(fid,'\t For build 1: [%3.1f %3.1f %3.1f %3.1f]\n',...
    countbest1(1),countbest1(2),countbest1(3),countbest1(4));
    fprintf(fid,'\t For build 2: [%3.1f %3.1f %3.1f %3.1f]\n',...
    countbest2(1),countbest2(2),countbest2(3),countbest2(4));
%
    countbest1m2 = (100/nkeepm2)*countbest1m2;
    countbest2m2 = (100/nkeepm2)*countbest2m2;
    fprintf(fid,'Without two-dimensional problems ');
    fprintf(fid,'(%d matrices):\n',nkeepm2);
    fprintf(fid,'\t For build 1: [%3.1f %3.1f %3.1f %3.1f]\n',...
    countbest1m2(1),countbest1m2(2),countbest1m2(3),countbest1m2(4));
    fprintf(fid,'\t For build 2: [%3.1f %3.1f %3.1f %3.1f]\n',...
    countbest2m2(1),countbest2m2(2),countbest2m2(3),countbest2m2(4));
%
    countbest1m3 = (100/nkeepm3)*countbest1m3;
    countbest2m3 = (100/nkeepm3)*countbest2m3;
    fprintf(fid,'Without two- and three-dimensional problems ');
    fprintf(fid,'(%d matrices):\n',nkeepm3);
    fprintf(fid,'\t For build 1: [%3.1f %3.1f %3.1f %3.1f]\n',...
    countbest1m3(1),countbest1m3(2),countbest1m3(3),countbest1m3(4));
    fprintf(fid,'\t For build 2: [%3.1f %3.1f %3.1f %3.1f]\n',...
    countbest2m3(1),countbest2m3(2),countbest2m3(3),countbest2m3(4));
end
%
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
