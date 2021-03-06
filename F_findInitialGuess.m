function initGuess = F_findInitialGuess(funct,materials,s,t,xmin,xmax,solParams)
%initGuess = F_findInitialGuess(funct,materials,s,t,xmin,xmax,solParams)
%computes the initial guess of the path for solving ultrasonic ray tracing
%with iterative root-finding algorithms.
%
%   Inputs:
%       funct - structured array of symbolic function and matlab functions,
%               as generated by the function: "F_genFunctions"
%       materials - structured array containing the proporties of the
%                   material of n layers
%           materials.v = [n x 1]  % wave propation speed in the material [mm/us]
%           materials.color = [n x 3]  % RGB color components for each material
%           materials.color = cell(n x 1)  % material names
%       s = [Sx,Sy] source coordinates
%       t = [Tx,Ty] target coordinates
%       xmin = a real number giving the left extremity of the interval
%              within which the intersection between the ray path and the
%              layer interfaces will be searched
%       xmax = a real number giving the right extremity of the interval
%              within which the intersection between the ray path and the
%              layer interfaces will be searched
%       solParameters = structured array containing all the parameters
%                       required for the definition of the algorithm to use
%           solParams.method = string containing the name of the algorithm
%                              to use. Possible values are: 'bisection',
%                              'bisectionWithNestedNR', 'NR', 'reducedNR',
%                              'hybridNR' and 'hybridReducedNR'.
%           solParams.numTol = Maximum numerical tolerance
%           solParams.stopMode = String indicating what stopping condition
%                                to use. Possible values are: 'dt', 'dx' 
%                                and 'dy'.
%           solParams.stopTol = Tolerance relative to the chosen stopping
%                               condition.
%   Outputs:
%       initGuess - structured array
%           initGuess.a = a real number giving the left extremity of the
%                         search interval for the bisection-based methods;
%           initGuess.b = a real number giving the right extremity of the
%                         search interval for the bisection-based methods;
%           initGuess.c = a real number giving the central value of the
%                         search interval for the bisection-based methods;
%           initGuess.fa = horizontal deviation in initGuess.a;
%           initGuess.fc = horizontal deviation in initGuess.c;
%
% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------

v = materials.v;
n = length(v);  %Number of layers
syms('x',[(2*n)-1,1]);

numTol = solParams.numTol;
method = solParams.method;

intF = funct.intF;
dIntF = funct.dIntF;
F = funct.F;
J = funct.J;

if strcmp(method,'bisectionWithNestedNR')
    newtonMode = true;
else
    newtonMode = false;
end


if n>1
    a0 = min(s(1),t(1)) - abs(t(1)-s(1))*0.5;
    b0 = max(s(1),t(1)) + abs(t(1)-s(1))*0.5;
    
    if (v(1)/v(2))<1
        fun = atan(dIntF(1)) - asin(v(1)/v(2)) - atan((x(1)-s(1))/(s(2)-intF(1)));
        fun = matlabFunction(fun,'Vars',x(1));
        a0 = F_bisectionMethod(fun,a0,b0,numTol);
        a = atan((a0-s(1))/(s(2)-double(subs(intF(1),'x1',a0))));
        
        fun = atan(dIntF(1)) + asin(v(1)/v(2)) - atan((x(1)-s(1))/(s(2)-intF(1)));
        fun = matlabFunction(fun,'Vars',x(1));
        b0 = F_bisectionMethod(fun,a0,b0,numTol);
        b = atan((b0-s(1))/(s(2)-double(subs(intF(1),'x1',b0))));
    else
        fun = atan(dIntF(1)) - (pi/2) - atan((x(1)-s(1))/(s(2)-intF(1)));
        fun = matlabFunction(fun,'Vars',x(1));
        a0 = F_bisectionMethod(fun,a0,b0,numTol);
        a = atan((a0-s(1))/(s(2)-double(subs(intF(1),'x1',a0))));
        
        fun = atan(dIntF(1)) + (pi/2) - atan((x(1)-s(1))/(s(2)-intF(1)));
        fun = matlabFunction(fun,'Vars',x(1));
        b0 = F_bisectionMethod(fun,a0,b0,numTol);
        b = atan((b0-s(1))/(s(2)-double(subs(intF(1),'x1',b0))));
    end
    
    if (abs(a-b)<=eps('single'))
        m = (t(2)-s(2))/(t(1)-s(1));
        c = s(2) - (m*s(1));
        fun = m*x(1) + c;
        fun = matlabFunction((fun - intF(1)),'Vars',x(1));
        a = min(s(1),t(1));
        b = max(s(1),t(1));
        root = F_bisectionMethod(fun,a,b,numTol);
        P1i = s(1) + ((v(1)/v(2))*(root-s(1)));
        
        if (v(1)/v(2))<1
            a = atan(double(subs(dIntF(1),'x1',P1i))) - asin(v(1)/v(2));
            b = atan(double(subs(dIntF(1),'x1',P1i))) + asin(v(1)/v(2));
        else
            a = atan(double(subs(dIntF(1),'x1',P1i))) - (pi/2);
            b = atan(double(subs(dIntF(1),'x1',P1i))) + (pi/2);
        end
    end
else
    a = atan((t(1)-s(1))/(s(2)-t(2)));
    b = a;
end

[~,fa] = F_propagateRay(intF,F,J,s,a,xmin,xmax,numTol,newtonMode);

if isnan(fa) && (a~=b)
    a0 = a;
    b0 = b;
    i = 0;
    while 1
        i = i+1;
        c = (a0+b0)/2;
        [~,fa] = F_propagateRay(intF,F,J,s,c,xmin,xmax,numTol,newtonMode);

        if isnan(fa)
            a0 = c;
        else
            b0 = c;
        end

        if i>10 && ~isnan(fa)
            a = c;
            break;
        end
        
        if i>100
            a = nan;
            b = nan;
            c = nan;
            fa = nan;
            fc = nan;
            break;
        end
    end
end

if ~isnan(a)
    [~,fb] = F_propagateRay(intF,F,J,s,b,xmin,xmax,numTol,newtonMode);
    
    if isnan(fb) && (a~=b)
        a0 = a;
        b0 = b;
        i = 0;
        while 1
            i = i+1;
            c = (a0+b0)/2;
            [~,fb] = F_propagateRay(intF,F,J,s,c,xmin,xmax,numTol,newtonMode);
            
            if isnan(fb)
                b0 = c;
            else
                a0 = c;
            end
            
            if i>10 && ~isnan(fb)
                b = c;
                break;
            end
        end
    end
    
    
    c = (a+b)/2;
    [c,fc] = F_propagateRay(intF,F,J,s,c,xmin,xmax,numTol,newtonMode);
end

initGuess.a = a;
initGuess.b = b;
initGuess.c = c;
initGuess.fa = fa;
initGuess.fc = fc;

%------------- END CODE --------------
end

