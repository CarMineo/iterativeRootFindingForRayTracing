function [propVars,deviation] = F_propagateRay(intF,F,J,s,iTheta,xmin,xmax,tol,newtonMode)
%[propVars,deviation] = F_propagateRay(intF,F,J,s,iTheta,xmin,xmax,tol,newtonMode)
%propagates the ray to all layers.
%
%   Inputs:
%       intF - symbolic interface functions
%       F = symbolic system of equations for the iterative methods;
%       J = symbolic Jacobian matrix relative to the system of equations;
%       s = [Sx,Sy] source coordinates
%       iTheta = input incidence angle with first interface
%       xmin = a real number giving the left extremity of the interval
%              within which the intersection between the ray path and the
%              layer interfaces will be searched
%       xmax = a real number giving the right extremity of the interval
%              within which the intersection between the ray path and the
%              layer interfaces will be searched
%       tol = Maximum numerical tolerance
%       newtonMode = boolean flag value to define if Newton-Raphson method
%                    is used (newtonMode = true) or if the bisection method
%                    is used (newtonMode = false) to find the intersection
%                    between the ray and the interfaces
%
%   Outputs:
%       propVar - propagation variables
%       deviation - horizontal deviation, measured on the horizontal line
%                   for the target.
%
% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------

n = (length(F)+1)/2;
syms('x',[(2*n)-1,1]);

propVars = nan(((2*n)-1),1);
propVars(n) = iTheta;

deviation = nan;

xVector = linspace(xmin,xmax,100);

for i=1:n-1
    intCurve = double(subs(intF(i),x(i),xVector));
    ymin = min(intCurve);
    ymax = max(intCurve);
    
    subSys = [F(i);F(n+i)];
    
    if i==1
        subSys = subs(subSys,x(n+i-1),propVars(n+i-1));
        if ymin~=ymax
            m = -1/tan(propVars(n+i-1));
            c = s(2) - (m*s(1));
            xmin = (ymax - c)/m;
            xmax = (ymin - c)/m;
        end
    else
        subSys = subs(subSys,{x(n+i-1);x(i-1)},{propVars(n+i-1);propVars(i-1)});
        if ymin~=ymax
            m = -1/tan(propVars(n+i-1));
            c = double(subs(intF(i-1),x(i-1),propVars(i-1))) - (m*propVars(i-1));
            xmin = (ymax - c)/m;
            xmax = (ymin - c)/m;
        end
    end
    
    fun = matlabFunction(subSys(1),'Vars',x(i));
    %[R,~,exitFlag] = fsolve(fun,0);
    
    if newtonMode == true
        xo = (xmin+xmax)/2;
        if i==1
            dFun = subs(J(i,i),x(n+i-1),propVars(n+i-1));
        else
            dFun = subs(J(i,i),{x(n+i-1);x(i-1)},{propVars(n+i-1);propVars(i-1)});
        end
        propVars(i) = F_newtonMethod(fun,dFun,xo,tol);
    else
        propVars(i) = F_bisectionMethod(fun,xmin,xmax,tol);
    end
    
    if ~isreal(propVars(i)) || isnan(propVars(i))
        break;
    end
    
    propVars(n+i) = double(subs(subSys(2),{x(i);x(n+i)},{propVars(i);0}));
end

if all(~isnan(propVars)) && (all(isreal(propVars)))
    if n==1
        deviation = double(subs(F(n),{x((2*n)-1)},{propVars((2*n)-1)}));
    else
        deviation = double(subs(F(n),{x((2*n)-1);x(n-1)},{propVars((2*n)-1);propVars(n-1)}));
    end
end

%------------- END CODE --------------

end

