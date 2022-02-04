function [minimum1,BestOrder1,minimum2,BestOrder2] = FindBestOrder(...
H,verboseflag,nocombi)
%% Test All Orders
if nocombi
%   Tries only 4 adhoc orders/heuristics
    P = [1;2;3;4];
else
%   Tries all possible orders and determine the best
    P = perms(1:size(H,1));
end
n=size(P,1);
%%
buildstyle=1;
for i = 1:n
    orderoption=P(i,:);
    [T,negativefound]=NES(H,orderoption,buildstyle,verboseflag);
    OrderResultBuild1(i)=negativefound;
end
minimum1=min(OrderResultBuild1);
BestOrder1=find(OrderResultBuild1==minimum1);
%fprintf('Build Type 1 results');
%eig(H)
%diag(H)
%P(BestOrder1,:)
%%
buildstyle=2;
for i = 1:n
    orderoption=P(i,:);
    [~,negativefound]=NES(H,orderoption,buildstyle,verboseflag);
    OrderResultBuild2(i)=negativefound;
end
minimum2=min(OrderResultBuild2);
BestOrder2=find(OrderResultBuild2==minimum2);
%fprintf('Build Type 2 results');
%eig(H);
%diag(H);
%P(BestOrder2,:);
%%
%if minimum1<minimum2
%    fprintf('\nThe overall best order was given by build type 1.\n')
%    fprintf('It occured %d different ways.\n',length(BestOrder1));
%    fprintf('It required %d updates.\n',minimum1);
%    fprintf('Whereas, build type 2 took %d updates.\n',minimum2);
%elseif minimum2<minimum1
%    fprintf('\nThe overall best order was given by build type 2.\n')
%    fprintf('It occured %d different ways.\n',length(BestOrder2));
%    fprintf('It required %d updates.\n',minimum2);
%    fprintf('Whereas, build type 1 took %d updates.\n',minimum1);
%else
%    fprintf('\nThe overall best order was given by build type 1 & 2.\n')
%    fprintf('It occured %d different type 1 ways.\n  ',length(BestOrder1));
%    fprintf('It occured %d different type 2 ways.\n',length(BestOrder2));
%    fprintf('It required %d updates.\n',minimum1);
%end
