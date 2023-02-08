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
%format long
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Hessian matrices
load('HESSIANS');
%
npbs = length(pbdims);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up
%
% Choose allowed dimensions 
maxdim = 50;
%maxdim = 16;
%maxdim=10;
%mindim=9;
%maxdim=8;
%maxdim = 5;
%maxdim = 11;
%mindim=10;
%mindim = 4;
mindim = 2;
%mindim = 16;
%mindim = 4;
% Select problems with appropriate dimensions
Imat = find(pbdims<=maxdim & pbdims>=mindim);
npbsel = length(pbdims(Imat)); % Number of problems selected
%
% Iterations to be kept
nitsNsel=nitsN;
% Finite difference values to be kept
%selFD = [1 2 3 4];
selFD = [1];
%selFD = [2 3 4];
%selFD = [4];
nFDsel = length(selFD);
% Total number of matrices
nmat = npbsel*(1+nitsNsel)*nFDsel;
%
% Select the ordering options (avoid combinatorial explosion)
if maxdim>8
    nocombi=1;
else
    nocombi=0;
end
%nocombi=0;
%nocombi=1;
%nocombi=-1;
% Plot all combinations in combinatorial case?
plotcombi=0;
%
% Random orthogonal transformation
randorthog=0;
%randorthog=1;
% Random permutation matrix
%randper=0;
randper=1;%
if randper && randorthog
    randorthog=0;
end
if randorthog || randper 
    rng(0);
end
%
verbose=0;%Verbose level in FindBestOrder subroutine
%
% Pre-allocate output data structures
%   minOrd* represents the minimal value of a certain ordering
%   bestOrd* contains the actual ordering that achieves the best value
minOrd1 = zeros(npbsel,nFDsel,1+nitsNsel);
bestOrd1 = cell(npbsel,nFDsel,1+nitsNsel);
minOrd2 = zeros(npbsel,nFDsel,1+nitsNsel);
bestOrd2 = cell(npbsel,nFDsel,1+nitsNsel);
%
%withnegdiag = 1;% Keep problems with negative diagonal elements?
withnegdiag = 0;% Remove problems with negative diagonal elements?
Ikeep = zeros(npbsel,nFDsel,1+nitsNsel);% 
%%%%%%%%%%%%%%%%%%%
% Main loop
for i=1:npbsel
%    auxi = nFDsel*(nitsNsel+1)*(i-1)+1;
    ipb=Imat(i);
    fprintf('Problem %s\n',pbnames{ipb});
    if randorthog
        myv = randn(pbdims(ipb),1);
        [myQ,~] = qr(myv);
    elseif randper
        myQ = eye(pbdims(ipb));
        myQ = myQ(randperm(pbdims(ipb)),:);
    end
%
%   Loop over finite difference values with the original matrix
    for iFD=1:nFDsel
        iselFD = selFD(iFD);
        valFD = hFD(iselFD);
        if valFD==0
            fprintf('\t Init pt (Exact) \n');
        else
            fprintf('\t Init Pt (FD=%1.0e) \n',valFD);
        end 
        if pbeigs(ipb,selFD(iFD),1)<-tolneg
            myH = pbmats{ipb}{iselFD}{1};
            if randorthog || randper
                myH = myQ*myH*myQ';
                if (sum(diag(myH)<0)>0)
                    mynegdiag=1;
                else
                    mynegdiag=0;
                end
            else
                mynegdiag=negdiags(ipb,iselFD,1);
            end 
            [minOrd1(i,iselFD,1),bestOrd1{i}{iselFD}{1},...
            minOrd2(i,iselFD,1),...
            bestOrd2{i}{iselFD}{1}] = FindBestOrder(myH,verbose,nocombi);
            if withnegdiag || ~mynegdiag
                Ikeep(i,iselFD,1) = 1;
            end
        else
            minOrd1(i,iselFD,1)=-1;
            bestOrd1{i}{iselFD}{1}=[];
            minOrd2(i,iselFD,1)=-1;
            bestOrd2{i}{iselFD}{1}=[];
        end
    end
    
%   Loop over the values with the Newton matrices
    for jN=1:nitsN
        for iFD=1:nFDsel
            iselFD = selFD(iFD);
            valFD = hFD(iselFD);
            if valFD==0
                fprintf('\t It Newton %d (Exact) \n',jN);
            else
                fprintf('\t It Newton %d (FD=%1.0e) \n',jN,valFD);
            end 
            if pbeigs(ipb,iselFD,1+jN)<-tolneg
                myH = pbmats{ipb}{iselFD}{1+jN};
                if randorthog || randper
                    myH = myQ*myH*myQ';
                    if (sum(diag(myH)<0)>0)
                        mynegdiag = 1;
                    else
                        mynegdiag = 0;
                    end
                else
                    mynegdiag = negdiags(ipb,iselFD,1+jN);
                end 
                [minOrd1(i,iselFD,1+jN),bestOrd1{i}{iselFD}{1+jN},...
                minOrd2(i,iselFD,1+jN),...
                bestOrd2{i}{iselFD}{1+jN}] = FindBestOrder(...
                myH,verbose,nocombi);
                if withnegdiag || ~mynegdiag
                    Ikeep(i,iselFD,1+jN) = 1;
                end
            else
                minOrd1(i,iselFD,1+jN)=-1;
                bestOrd1{i}{iselFD}{1+jN}=[];
                minOrd2(i,iselFD,1+jN)=-1;
                bestOrd2{i}{iselFD}{1+jN}=[];
            end
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%
%
% Save the relevant data into a .mat file
save DATANESCUTEST pbdims nitsN selFD Imat Ikeep minOrd1 bestOrd1 minOrd2 bestOrd2
%
%%%%%%%%%%%%%%%%%%%
% Write the desired outputs in a data file
fid = fopen('ResultsHessianCUTEst','w');
fprintf(fid,'Minimum allowed dimension: %d\n\n',mindim);
fprintf(fid,'Maximum allowed dimension: %d\n\n',maxdim);
switch nocombi
    case 1
        fprintf(fid,'Orderings:\n');
        fprintf(fid,'\t 1 - 1:n\n');
        fprintf(fid,'\t 2 - Smallest diag to largest diag\n');
        fprintf(fid,'\t 3 - Largest diag to smallest diag\n');
        fprintf(fid,'\t 4 - Flipflop smallest/largest diag\n\n');
    case -1
        fprintf(fid,'Ordering:\n');
        fprintf(fid,'\t 1 - 1:n\n');
    otherwise
    fprintf(fid,'Orderings: Combinatorial\n\n');
end
% Determine winners between Build 1 and Build 2
%nkeep = length(Ikeep);
auxkeep = find(Ikeep==1);
%size(auxkeep)
nkeep = sum(Ikeep,'all');
%selFD;
nkeepm2 = 0;
nkeepm3 = 0;
bestmin1 = (minOrd1<minOrd2);
sb1 = sum(bestmin1(auxkeep),'all');
bestmin2 = (minOrd1>minOrd2);
sb2 = sum(bestmin2(auxkeep),'all');
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
switch nocombi
    case 1
        countbest1=[0 0 0 0];
        countbest1m2 = [0 0 0 0];
        countbest1m3 = [0 0 0 0];
        countbest2=[0 0 0 0];
        countbest2m2 = [0 0 0 0];
        countbest2m3 = [0 0 0 0];
    case -1
        countbest1=[0];
        countbest1m2 = [0];
        countbest1m3 = [0];
        countbest2=[0];
        countbest2m2 = [0];
        countbest2m3 = [0];
end
%
fprintf(fid,'Problem \t| Dim | It Newton | FinDiff | MinEig |');
fprintf(fid,' BestBuild | MinUpdate | BestOrder\n');
fprintf(fid,'-----------------------------------------------------\n\n');
%%%%%%%%%
% Detailed results for all problems
% Outer loop on CUTEst problems, inner loop on the various matrices 
% with respect to that problem.
for i=1:npbsel
    ipb=Imat(i);
    switch nocombi
        case 0
            if plotcombi
                Pi = perms(1:pbdims(ipb));
            end
            nways=factorial(pbdims(ipb));
        case 1
            nways=4;
        case -1
            nways=1;
    end
%    auxi = (nFD+1)*(nitsN+1)*(i-1);
%    auxj = auxi + 1;
    for j=0:1:nitsN
        for iFD=1:nFDsel
            jFD = selFD(iFD);
            fprintf('Problem %s ItN %d FD %1.0e\n',...
            pbnames{ipb},j,hFD(jFD));
            fprintf(fid,'%s \t %d \t %d \t',pbnames{ipb},...
            pbdims(ipb),j);
            if jFD==0
                fprintf(fid,'Exact \t');
            else
                fprintf(fid,'%1.0e \t',hFD(jFD));
            end
            fprintf(fid,'%1.3e \t\t',pbeigs(ipb,jFD,1+j));
            if Ikeep(i,jFD,1+j)
                if pbdims(ipb)>=3
                    nkeepm2 = nkeepm2+1;
                    if pbdims(ipb)>3
                        nkeepm3 = nkeepm3+1;
                    end
                end
                if bestmin1(i,jFD,1+j)
                    fprintf(fid,'1 \t %d \t',minOrd1(i,jFD,1+j));
                    for k=1:length(bestOrd1{i}{jFD}{1+j})
                        if (nocombi~=0) || (~nocombi & plotcombi)
                            fprintf(fid,'%d, ',bestOrd1{i}{jFD}{1+j}(k));
                            if (nocombi~=0)
                                countbest1(...
                                bestOrd1{i}{jFD}{1+j}(k)) = countbest1(...
                                bestOrd1{i}{jFD}{1+j}(k))+1;
                                if pbdims(ipb)>=3
                                    countbest1m2(...
                                    bestOrd1{i}{jFD}{1+j}(k)) = countbest1m2(...
                                    bestOrd1{i}{jFD}{1+j}(k))+1;
                                    if pbdims(ipb)>3
                                        countbest1m3(...
                                        bestOrd1{i}{jFD}{1+j}(k)) = countbest1m3(...
                                        bestOrd1{i}{jFD}{1+j}(k))+1;
                                    end
                                end
                            end
                        else
                            fprintf(fid,'[ ');
                            for l=1:pbdims(ipb)
                                fprintf(fid,'%d ',Pi(...
                                bestOrd1{i}{jFD}{1+j}(k),l));
                            end
                            fprintf(fid,'],');
                        end
                    end                
                    fprintf(fid,' \n');
                else
                    if bestmin2(i,jFD,1+j)
                        fprintf(fid,'2 \t %d \t',minOrd2(i,jFD,1+j));
                        for k=1:length(bestOrd2{i}{jFD}{1+j})
                            if (nocombi~=0) || (~nocombi & plotcombi)
                                fprintf(fid,'%d, ',...
                                bestOrd2{i}{jFD}{1+j}(k));
                                if (nocombi~=0)
                                    countbest2(...
                                    bestOrd2{i}{jFD}{1+j}(k)) = countbest2(...
                                    bestOrd2{i}{jFD}{1+j}(k))+1;
                                    if pbdims(ipb)>=3
                                        countbest2m2(...
                                        bestOrd2{i}{jFD}{1+j}(k)) = countbest2m2(...
                                        bestOrd2{i}{jFD}{1+j}(k))+1;
                                        if pbdims(ipb)>3
                                            countbest2m3(...
                                            bestOrd2{i}{jFD}{1+j}(k)) = countbest2m3(...
                                            bestOrd2{i}{jFD}{1+j}(k))+1;
                                        end
                                    end
                                end
                            end
                        end
                        fprintf(fid,' \n');
                    else
                        fprintf(fid,'1+2 \t %d \t',minOrd1(i,jFD,1+j));
                        for k=1:length(bestOrd1{i}{jFD}{1+j})
                            if (nocombi~=0) || (~nocombi & plotcombi)
                                fprintf(fid,'%d, ',bestOrd1{i}{jFD}{1+j}(k));
                                if (nocombi~=0)
                                    countbest1(...
                                    bestOrd1{i}{jFD}{1+j}(k)) = countbest1(...
                                    bestOrd1{i}{jFD}{1+j}(k))+1;
                                    if pbdims(ipb)>=3
                                        countbest1m2(...
                                        bestOrd1{i}{jFD}{1+j}(k)) = countbest1m2(...
                                        bestOrd1{i}{jFD}{1+j}(k))+1;
                                        if pbdims(ipb)>3
                                            countbest1m3(...
                                            bestOrd1{i}{jFD}{1+j}(k)) = countbest1m3(...
                                            bestOrd1{i}{jFD}{1+j}(k))+1;
                                        end
                                    end
                                end
                            end
                        end
                        fprintf(fid,'+');
                        for k=1:length(bestOrd2{i}{jFD}{1+j})
                            if (nocombi~=0) || (~nocombi & plotcombi)
                                fprintf(fid,'%d, ',bestOrd2{i}{jFD}{1+j}(k));
                                if (nocombi~=0)
                                    countbest2(...
                                    bestOrd2{i}{jFD}{1+j}(k)) = countbest2(...
                                    bestOrd2{i}{jFD}{1+j}(k))+1;
                                    if pbdims(ipb)>=3
                                        countbest2m2(...
                                        bestOrd2{i}{jFD}{1+j}(k)) = countbest2m2(...
                                        bestOrd2{i}{jFD}{1+j}(k))+1;
                                        if pbdims(ipb)>3
                                            countbest2m3(...
                                            bestOrd2{i}{jFD}{1+j}(k)) = countbest2m3(...
                                            bestOrd2{i}{jFD}{1+j}(k))+1;
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
    end
end
%
%%%%%%%%%%%%%%%%%%%
% Final statistics:
% For dimension-independent strategies, plot the percentage of pbms for 
% which each strategy was the best
fprintf(fid,'=============================================================\n');
if (nocombi==1)
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
