sigma = 5;
b = 3;
N = 10000;
[c,u] = myQuantizer(sigma, b, N);% centres and borders


function [c,u] = myQuantizer(sigma, b, N)
    X = normrnd(0,sigma,[1,N]); % N random numbers with normal distribution
    X = sort(X);
    u = datasample(X, (2.^b)-1, 'Replace', false); % initialize random borders
    u = sort(u);
    u_prev = zeros(size(u));
    [c,u] = update_borders(X, u, b, N);
    iters = 0;
    while  max(abs(u-u_prev))>0.0002 % check convergence
        iters = iters + 1;
        u_prev = u;
        [c,u] = update_borders(X, u, b, N);
    end
    
end

function [newC, newU] = update_borders(X, u, b, N)
    newC = zeros(2.^b, 1);
    newU = zeros(1, (2.^b)-1);
    for n = 1:2^b
         if n==1
              [row, col] = find(X<=u(1,1));
         elseif n==(2^b)
              [row, col] = find(X>u(1,(2.^b)-1));
         else
              [row, col] = find((u(1,n-1)<X) & (X<=u(1,n)));
        end
        G=X(1,col);       
        norms = zeros(1,length(col));
        for j = 1:length(col)
            norms(1,j) = sum((repmat(G,1,1) - X(1,col(j))) .^ 2, 2);
        end
        [M,I] = min(norms);
        newC(n,1)= X(1,col(I));
    end
    for n = 1:(2^b)-1
        newU(1,n) = (newC(n,1)+newC(n+1,1))/2;
    end
    
end
