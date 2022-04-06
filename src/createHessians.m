%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting Hessian matrices from CUTest problems
%
% This script calls the CUTEst/MATLAB interface to obtain Hessian at initial 
% points of CUTEst problems, and finite-differences approximations thereof.
%
% In the same directory, it is assumed that there exists a file ListPbmsNC 
% containing the data related to every problem in the collection under the 
% form:
%
% NameProblem	Dimension   Best Known value (not used for this script)
%
% An data output file is produced that contains the following cell structures:
%   pbnames: cell structure that contains the names of the CUTEst 
%   problems that were used.
%   pbdims: vector structure containing the corresponding dimensions.
%   pbmats: cell structure containing matrices. These matrices are 
%   the Hessian matrices computed at the initial points of every CUTEst 
%   problem in the test set.
%   pbeigs: vector structure containing the minimum eigenvalues of every 
%   matrix in pbmats (these eigenvalues are computed via eigs).
%   negdiags: vector structure of same length than pbeigs that contains 
%   0/1 values, indicating whether the corresponding matrix in pbmats has a 
%   negative diagonal element.
%
% Other optional fields may be saved:
%   pbmatsN: cell structure containing matrices. These matrices are 
%   the Hessian matrices computed at points obtaining by applying Newton's 
%   method to  initial points of every CUTEst 
%   problem in the test set.
%   pbeigsN: matrix structure containing the minimum eigenvalues of every 
%   matrix in pbmatsN (these eigenvalues are computed via eigs).
%   negdiagsN: vector structure of same length than pbeigs that contains 
%   0/1 values, indicating whether the corresponding matrix in pbmats has a 
%   negative diagonal element.
%   pbmatsFD: cell structure containing matrices. These matrices are 
%   the Hessian matrices computed at the initial points of every CUTEst 
%   problem in the test set.
%   pbeigsFD: vector structure containing the minimum eigenvalues of every 
%   matrix in pbmats (these eigenvalues are computed via eigs).
%   negdiagsFD: vector structure of same length than pbeigs that contains 
%   0/1 values, indicating whether the corresponding matrix in pbmats has a 
%   negative diagonal element.
%   pbmatsNFD: cell structure containing matrices. These matrices are 
%   the Hessian matrices computed at the initial points of every CUTEst 
%   problem in the test set.
%   pbeigsNFD: vector structure containing the minimum eigenvalues of every 
%   matrix in pbmats (these eigenvalues are computed via eigs).
%   negdiagsNFD: vector structure of same length than pbeigs that contains 
%   0/1 values, indicating whether the corresponding matrix in pbmats has a 
%   negative diagonal element.
%
% Implementation: C. W. Royer 
% Started July 29, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%
% Secure path to the CUTEst/MATLAB interface
addpath('/home/croyer/Documents/Codes/cutest/cutest/src/matlab')
% Getting the problem list
fid = fopen('ListPbmsNC','r');
%fid = fopen('ListOnePb','r');
PROBS = textscan(fid,'%s %d %f');
fclose(fid);
%fid2 = fopen('HessOnePb','w');
%
npbs=length(PROBS{1});
pbnames = PROBS{1};
pbdims = PROBS{2};
BVALS = PROBS{3};
%
findiff=1;%Using finite difference approximation to the Hessian matrix
hFD = [1e-2,1e-4,1e-6];
nh = length(hFD);
findiff = (nh>0);
% Creating output file and structure
%
if findiff
    fid2 = fopen('HessPbmsNCFinDiff','w');
else
    fid2 = fopen('HessPbmsNC','w');
end
fprintf(fid2,'Problem Name & Dimension & InitPoint (Exact) ');
if findiff
    for i=1:nh
        fprintf(fid2,' & InitPoint (FD=%1.2e)',hFD(i));
    end
end
% Defining structures to store matrices and eigenvalue information
pbmats = cell(npbs,1);
pbeigs = zeros(npbs,1);
negdiags = zeros(npbs,1);
% Optional: Use finite-differences to approximate the matrices
if findiff
    pbmatsFD = cell(npbs,nh);
    pbeigsFD = zeros(npbs,nh);
    negdiagsFD = zeros(npbs,nh);
end
%
% Optional: Perform Newton's iteration to get more Hessian matrices
nitsN = 2;
for i=1:nitsN
    fprintf(fid2,'& It %d Newton (Exact)',i);
%   Optional: Use finite differences to approximate those matrices
    if findiff
        for j=1:nh
            fprintf(fid2,'& It %d Newton (FD=%1.2e)',i,hFD(j));
        end
    end
end
fprintf(fid2,'\n\n');
if nitsN>0
    pbmatsN = cell(npbs,nitsN);
    pbeigsN = zeros(npbs,nitsN);
    negdiagsN = zeros(npbs,nitsN);
    if findiff
        pbmatsNFD = cell(npbs,nitsN,nh);
        pbeigsNFD = zeros(npbs,nitsN,nh);
        negdiagsNFD = zeros(npbs,nitsN,nh);
    end
else
    pbmatsN = [];
    pbeigsN = [];
    if findiff
        pbmatsNFD = [];
        pbeigsNFD = [];
    end
end
%
% Number of matrices with negative diagonal elements
nb_negdiag=0;
nb_negcurv=0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for numpb = 1:npbs
%
%	Decoding the SIF file of the problem
   	st=(['!cutest2matlab ' PROBS{1}{numpb}]);
   	disp(st); eval(st); pause(1);
%
%	Gathering the problem data
   	prob = cutest_setup();
   	x0=prob.x;
	n= prob.n;
   	name=prob.name;
%
	if (pbdims(numpb)~=n)
		error('CUTEST:pb %s has wrong dimension',name);
	end
%
	fprintf('%s & %d\n',name,n);
    g0 = cutest_grad(x0);
	H0 = cutest_hess(x0);
    pbmats{numpb}=H0;
    if (sum(diag(H0)<0)>0)
        negdiags(numpb)=1;
        nb_negdiag=nb_negdiag+1;
    end
	[~,l0] = eigs(H0,1,'sa');
    pbeigs(numpb)=l0;
    if l0<0
        nb_negcurv=nb_negcurv+1;
    end
	fprintf(fid2,'%s & %d & %1.3e (%d) ',name,n,l0,negdiags(numpb));
    if findiff
        for iFD=1:nh
            [H0FD,flagFD] = findiffapprox(x0,@(x) cutest_obj(x),hFD(iFD));
            pbmatsFD{numpb}{iFD}=H0FD;
            if ~flagFD
                if (sum(diag(H0FD<0))>0)
                    negdiagsFD(numpb,iFD)=1;
                    nb_negdiag=nb_negdiag+1;
                end
                [~,l0FD] = eigs(H0FD,1,'sa');
                pbeigsFD(numpb,iFD)=l0FD;
                if l0FD<0
                    nb_negcurv=nb_negcurv+1;
                end
            else
                pbeigsFD(numpb,iFD)=NaN; 
            end
	        fprintf(fid2,'& %1.3e (%d) ',l0FD,negdiagsFD(numpb,iFD));
        end
    end
%
%   Optional: perform iterations of Newton's method to collect  more matrices
%
    for i=1:nitsN
        x0 = x0 - H0 \ g0;
        g0 = cutest_grad(x0);
        H0 = cutest_hess(x0);
        pbmatsN{numpb}{i}=H0;
        if (sum(diag(H0)<0)>0)
            negdiagsN(numpb,i)=1;
            nb_negdiag=nb_negdiag+1;
        end
	    [~,l0] = eigs(H0,1,'sa');
        if l0<0
            nb_negcurv=nb_negcurv+1;
        end
        pbeigsN(numpb,i)=l0;
	    fprintf(fid2,'& %1.3e (%d) ',l0,negdiagsN(numpb,i));
        if findiff
            for iFD=1:nh
                [H0FD,flagFD] = findiffapprox(x0,@(x) cutest_obj(x),hFD(iFD));
                pbmatsNFD{numpb}{i}{iFD}=H0FD;
                if ~flagFD
                    if (sum(diag(H0FD)<0)>0)
                        negdiagsNFD(numpb,i,iFD)=1;
                        nb_negdiag=nb_negdiag+1;
                    end
                    [~,l0FD] = eigs(H0FD,1,'sa');
                    pbeigsNFD(numpb,i,iFD)=l0FD;
                    if l0FD<0
                        nb_negcurv=nb_negcurv+1;
                    end
                else
                    pbeigsNFD(numpb,i,iFD)=NaN;
                end
	            fprintf(fid2,'& %1.3e (%d) ',l0FD,negdiagsNFD(numpb,i,iFD));
            end
        end
    end
    fprintf(fid2,'\n');
%
end
%
nbmats_total = (nitsN+1)*(nh+1)*npbs;
fprintf('Matrices with negative diagonal elements: %d out of %d\n',...
nb_negdiag,nbmats_total);
fprintf('Matrices with negative curvature: %d out of %d\n',...
nb_negcurv,nbmats_total);
%
fclose(fid2);
%
if ~findiff
    save HESSIANSNC pbnames pbdims pbmats pbeigs nitsN pbmatsN pbeigsN findiff negdiags negdiagsN nb_negdiag nb_negcurv
else
    save HESSIANSNCFD pbnames pbdims pbmats pbeigs nitsN pbmatsN pbeigsN hFD pbmatsFD pbeigsFD pbmatsNFD pbeigsNFD findiff negdiags negdiagsFD negdiagsN negdiagsNFD nb_negdiag nb_negcurv
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
