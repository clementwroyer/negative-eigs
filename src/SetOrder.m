function Order = SetOrder(D,n,orderoption,buildstyle)
%% SetOrder, by Warren Hare, 2021
% Used in NES to create a 2 x (n^2-n)/2 matrix that orders the growth of
% Hestimate
% Orders: vector = use the order given by the vector
%         1 = 1, 2, 3, ..., n
%         2 = smallest diag entry to largest diag entry
%         3 = largest diag entry to smallest diag entry
%         4 = flip-flop smallest, largest, 2nd smallest, 2nd largest, etc
% Buildstyle: 1 = scatter the matrix building many small submatrices
%             2 = focuses on building up large submatrices
%
    if length(orderoption)>1
        I = orderoption;
    elseif orderoption==1
        I=1:n; % This order is 1, 2, ..., n
    elseif orderoption==2
        [~,I]=sort(D);I=I'; % This order is the smallest to largest from the diagonal
    elseif orderoption==3
        [~,I]=sort(D,'descend');I=I'; % This order is the largest to smallest from the diagonal
    elseif orderoption==4
        [~,temp1]=sort(D);
        [~,temp2]=sort(D,'descend');
        if (sum(temp1==temp2)==n) || (n<=2)
            I = temp1;
        else
            for i=1:n
                if mod(i,2) 
                    I(i)=temp1((i+1)/2);
                else
                    I(i)=temp2(i/2);
                end
            end
        end
    end
    Order=[];
    I;
    % This scatter the matrix building many small submatrices
    if buildstyle==1
        for counti=1:n
            for countj=counti+1:n
                Order=[Order,[I(counti);I(countj)]];
            end
        end
    else
    % This focuses on building up large submatrices
        for counti=2:n
            for countj=counti-1:-1:1
                Order=[Order,[I(counti);I(countj)]];
            end
        end
    end
    %Order
end
