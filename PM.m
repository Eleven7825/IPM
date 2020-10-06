%power method

function [mu,x,k] = PM(A,x,TOL,N)
%input parameters

%set up iteration time k and normalize x
k = 1;
x = x/norm(x);

%Do the iteration
while k < N + 1
    y = A*x;
    mu = norm(y);
    if y == 0
        disp('Eigenvector:')
        disp(x)
        disp('A has the eigenvalue 0, select a new vector x and restart')
        x = input('input the initial vector:');
    end
    y = y / norm(y);
    
    ERR = abs(norm(A*y) - mu);
    if ERR < TOL 
        return
    end
    x = y;
    k = k + 1;
    if k == N
      error('Max iteration time exceeded! ')
    end
end
end

