function [root] = F_bisectionMethod(fun,a,b,tol)
%[root] = F_bisectionMethod(fun,a,b,tol) performs the bisection algorithm
%to find the root of the function given by "fun", within the interval [a,b].
%
%   Inputs:
%       fun - single variable function
%       a - a real number giving the left extremity of the search interval
%       b - a real number giving the right extremity of the search interval
%       tol - tolerance to stop the iteration of the computation
%
%   Outputs:
%       root - root of the function in the interval [a,b]
%
% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------

co = 0;
cn = inf;
fa = fun(a);

while abs(cn-co) > tol
    co = cn;
    cn = (a+b)/2;
    fc = fun(cn);
    if sign(fc)~=sign(fa)
        b = cn;
    else
        if abs(fc)<abs(fa)
            a = cn;
            fa = fc;
        else
            b = cn;
        end 
    end
end

root = cn;

%------------- END CODE --------------

end

