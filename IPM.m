% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of 
% the GNU Lesser General Public License, either version 3 of the License, 
% or any later version.

function [eig_value,eig_vector,i] = IPM(A,mu,x0,m,e)
% initialize the method 
x = x0;
B = A - mu*eye(size(A));
i = 0;
ERR = 1+e;
Xp = 0;
for j = 1:length(x)
    if abs(x(j)) > Xp
        Xp = x(j);
    end
end
x = x/Xp;

while i < m && ERR > e
    y = B\x;
    [~,loc] = max(abs(y));
    Yp = y(loc);
    
    ERR = 0;
    if Yp == 0
        Yp = -1/mu;
        break
    end
    
    errvec = x - (y/Yp);
    for j = 1:length(errvec)
        if abs(errvec(j)) > ERR
            ERR = abs(errvec(j));
        end
    end

    x = y/Yp;
    i = i+1;
end

if i < m
    eig_value = 1/Yp+mu;
    eig_vector = x;
    return
else
    eig_value = NaN;
    eig_vector = NaN;
    return
end
end
