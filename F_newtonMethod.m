function [root] = F_newtonMethod(fun,dFun,xo,tol)
%[root] = F_newtonMethod(fun,dFun,xo,tol) performs the Newton-Raphson algorithm
%to find the root of the function given by "fun", using xo as initial guess.
%
%   Inputs:
%       fun - single variable function
%       dFun - derivative of single variable function
%       xo - initial guess
%       tol - tolerance to stop the iteration of the computation
%
%   Outputs:
%       root - root of the function
%
% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------

dFun = matlabFunction(dFun);
xn = xo - (fun(xo)/dFun(xo));
while abs(xn-xo) > tol
    xo = xn;
    xn = xo - (fun(xo)/dFun(xo));
end
root = xn;

%------------- END CODE --------------

end

