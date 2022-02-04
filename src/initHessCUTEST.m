%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting Hessian matrices from CUTest problems
%
% This script calls the CUTEst/MATLAB interface to obtain Hessian at initial 
% points of CUTEst problems.
%
% In the same directory, it is assumed that there exists a file ListPbmsNC 
% containing the data related to every problem in the collection under the 
% form:
%
% NameProblem	Dimension   Best Known value (not used for this script)
%
% An data output file is produced that contains the following cell structures:
%   pbnames is a cell structure that contains the names of the CUTEst 
%   problems that were used.
%   pbdims is a vector structure containing the corresponding dimensions.
%   pbmats is a cell structure containing matrices. These matrices are 
%   the Hessian matrices computed at the initial points of every CUTEst 
%   problem in the test set.
%   pbeigs is a vector structure containing the minimum eigenvalues of every 
%   matrix in pbmats (these eigenvalues are computed via eigs).
%
%
% Implementation: C. W. Royer 
% Started July 29, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%
% Secure path to the CUTEst/MATLAB interface
% Swann
%addpath('/home/croyer/Documents/TheCUTEst/cutest/src/matlab')
% Rastignac
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
pbmats = cell(npbs,1);
pbeigs = zeros(npbs,1);
if findiff
    pbmatsFD = cell(npbs,nh);
    pbeigsFD = zeros(npbs,nh);
end
%
% Optional: Perform Newton's iteration to get more Hessian matrices
nitsN = 2;
for i=1:nitsN
    fprintf(fid2,'& It %d Newton (Exact)',i);
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
    if findiff
        pbmatsNFD = cell(npbs,nitsN,nh);
        pbeigsNFD = zeros(npbs,nitsN,nh);
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
	[~,l0] = eigs(H0,1,'sa');
    pbeigs(numpb)=l0;
	fprintf(fid2,'%s & %d & %1.3e ',name,n,l0);
    if findiff
        for iFD=1:nh
            [H0FD,flagFD] = findiffapprox(x0,@(x) cutest_obj(x),hFD(iFD));
            pbmatsFD{numpb}{iFD}=H0FD;
            if ~flagFD
                [~,l0FD] = eigs(H0FD,1,'sa');
                pbeigsFD(numpb,iFD)=l0FD;
            else
                pbeigsFD(numpb,iFD)=NaN; 
            end
	        fprintf(fid2,'& %1.3e ',l0FD);
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
	    [~,l0] = eigs(H0,1,'sa');
        pbeigsN(numpb,i)=l0;
	    fprintf(fid2,'& %1.3e ',l0);
        if findiff
            for iFD=1:nh
                [H0FD,flagFD] = findiffapprox(x0,@(x) cutest_obj(x),hFD(iFD));
                pbmatsNFD{numpb}{i}{iFD}=H0FD;
                if ~flagFD
                    [~,l0FD] = eigs(H0FD,1,'sa');
                    pbeigsNFD(numpb,i,iFD)=l0FD;
                else
                    pbeigsNFD(numpb,i,iFD)=NaN;
                end
	            fprintf(fid2,'& %1.3e ',l0FD);
            end
        end
    end
    fprintf(fid2,'\n');
%
end
%
fclose(fid2);
%
if ~findiff
    save HESSIANSNC pbnames pbdims pbmats pbeigs nitsN pbmatsN pbeigsN findiff
else
    save HESSIANSNCFD pbnames pbdims pbmats pbeigs nitsN pbmatsN pbeigsN hFD pbmatsFD pbeigsFD pbmatsNFD pbeigsNFD findiff
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
