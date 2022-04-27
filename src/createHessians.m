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
% An data output file HESSIANS is produced that contains the following cell 
% structures:
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
%
npbs=length(PROBS{1});
pbnames = PROBS{1};
pbdims = PROBS{2};
BVALS = PROBS{3};
%
%Using finite difference approximation to the Hessian matrix
hFD = [0,1e-2,1e-4,1e-6];
%hFD=[0];
%hFD = [];
nh = length(hFD);
% Optional: Build more matrices by performing iterations of Newton's method
nitsN = 2;
% Creating output file and structure
%
fid2 = fopen('HessianEigs','w');
fprintf(fid2,'Problem Name & Dimension ');
for i=1:nh
    if hFD(i)==0
        fprintf(fid2,' & InitPoint (Exact) ');
    else
        fprintf(fid2,' & InitPoint (FD=%1.2e)',hFD(i));
    end
end
for j=1:nitsN
    for i=1:nh
        if hFD(i)==0
            fprintf(fid2,' & It %d Newton (Exact) ',j);
        else
            fprintf(fid2,' & It %d Newton (FD=%1.2e)',j,hFD(i));
        end
    end
end
fprintf(fid2,'\n\n');
% Defining structures to store matrices and eigenvalue information
pbmats = cell(npbs,nh,1+nitsN);
pbeigs = zeros(npbs,nh,1+nitsN);
negdiags = zeros(npbs,nh,1+nitsN);
fprintf(fid2,'\n\n');
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

%   Compute matrices related to the initial point

    for iFD=1:nh
        if hFD(iFD)==0
            pbmats{numpb}{iFD}{1}=H0;
            if (sum(diag(H0)<0)>0)
                negdiags(numpb,iFD,1)=1;
                nb_negdiag=nb_negdiag+1;
            end
	        [~,l0] = eigs(H0,1,'sa');
            pbeigs(numpb,iFD,1)=l0;
            if l0<0
                nb_negcurv=nb_negcurv+1;
            end
            fprintf(fid2,'%s & %d & %1.3e (%d) ',name,n,l0,...
            negdiags(numpb,iFD,1));
        else
            [H0FD,flagFD] = findiffapprox(x0,@(x) cutest_obj(x),hFD(iFD));
            pbmats{numpb}{iFD}{1}=H0FD;
            if ~flagFD
                if (sum(diag(H0FD<0))>0)
                    negdiags(numpb,iFD,1)=1;
                    nb_negdiag=nb_negdiag+1;
                end
                [~,l0FD] = eigs(H0FD,1,'sa');
                pbeigs(numpb,iFD,1)=l0FD;
                if l0FD<0
                    nb_negcurv=nb_negcurv+1;
                end
            else
                pbeigs(numpb,iFD,1)=NaN; 
            end
	        fprintf(fid2,'& %1.3e (%d) ',l0FD,...
            negdiags(numpb,iFD,1));
        end
    end
%
%   Optional: perform iterations of Newton's method to collect more matrices
%
    for jN=1:nitsN
        x0 = x0 - H0 \ g0;
        g0 = cutest_grad(x0);
        H0 = cutest_hess(x0);

%             
        for iFD=1:nh

            if hFD(iFD)==0
%               Exact matrix
                pbmats{numpb}{iFD}{1+jN}=H0;
                if (sum(diag(H0)<0)>0)
                    negdiags(numpb,iFD,1+jN)=1;
                    nb_negdiag=nb_negdiag+1;
                end
	            [~,l0] = eigs(H0,1,'sa');
                if l0<0
                    nb_negcurv=nb_negcurv+1;
                end
                pbeigs(numpb,iFD,1+jN)=l0;
                fprintf(fid2,'& %1.3e (%d) ',l0,negdiags(numpb,iFD,1+jN));
            else
%               Finite-difference formula
                [H0FD,flagFD] = findiffapprox(x0,@(x) cutest_obj(x),hFD(iFD));
                pbmats{numpb}{iFD}{1+jN}=H0FD;
                if ~flagFD
                    if (sum(diag(H0FD)<0)>0)
                        negdiags(numpb,iFD,1+jN)=1;
                        nb_negdiag=nb_negdiag+1;
                    end
                    [~,l0FD] = eigs(H0FD,1,'sa');
                    pbeigs(numpb,iFD,1+jN)=l0FD;
                    if l0FD<0
                        nb_negcurv=nb_negcurv+1;
                    end
                else
                    pbeigs(numpb,iFD,1+jN)=NaN;
                end
	            fprintf(fid2,'& %1.3e (%d) ',l0FD,...
                negdiags(numpb,iFD,1+jN));
            end
        end
    end
    fprintf(fid2,'\n');
%
end
%
nbmats_total = (nitsN+1)*nh*npbs;
fprintf('Matrices with negative diagonal elements: %d out of %d\n',...
nb_negdiag,nbmats_total);
fprintf('Matrices with negative curvature: %d out of %d\n',...
nb_negcurv,nbmats_total);
%
fclose(fid2);
%
save HESSIANS pbnames pbdims pbmats pbeigs nitsN hFD negdiags nb_negdiag nb_negcurv
%if ~findiff
%   save HESSIANS pbnames pbdims pbmats pbeigs nitsN pbmatsN pbeigsN findiff negdiags negdiagsN nb_negdiag nb_negcurv
%else
%    save HESSIANS pbnames pbdims pbmats pbeigs nitsN pbmatsN pbeigsN hFD pbmatsFD pbeigsFD pbmatsNFD pbeigsNFD findiff negdiags negdiagsFD negdiagsN negdiagsNFD nb_negdiag nb_negcurv
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
