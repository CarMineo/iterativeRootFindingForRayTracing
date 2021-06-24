function funct = F_genFunctions(intF,materials,source,target,solParams)
%funct = F_genFunctions(intF,materials,source,target,solParams)
%generates the symbolic functions for the iterative root-finding algorithms
%
%   Inputs:
%       intF - array of n-1 symbolic functions of the interfaces, for n
%              layers
%       materials - structured array containing the proporties of the
%                   material of n layers
%           materials.v = [n x 1]  % wave propation speed in the material [mm/us]
%           materials.color = [n x 3]  % RGB color components for each material
%           materials.color = cell(n x 1)  % material names
%       source = [Sx,Sy] source coordinates
%       target = [Tx,Ty] target coordinates
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
%       funct - structured array
%           funct.intF = a copy of the array containing the symbolic 
%                        interface functions;
%           funct.dIntF = symbolic derivatives of the interface functions;
%           funct.F = symbolic system of equations for the iterative methods;
%           funct.J = symbolic Jacobian matrix relative to the system of equations;
%           funct.FmF = symbolic matlab functions of the system of equations for the iterative methods;
%           funct.JmF = symbolic matlab functions of the Jacobian matrix;
%           funct.rF = symbolic reduced system of equations for the iterative methods;
%           funct.rJ = symbolic Jacobian matrix relative to the reduced system of equations;
%           funct.rFmF = symbolic matlab functions of the reduced system of equations for the iterative methods;
%           funct.rJmF = symbolic matlab functions of the Jacobian matrix relative to the reduced system;
%
% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------

funct = [];

funct.intF = intF;
funct.dIntF = [];
funct.F = [];
funct.J = [];
funct.FmF = [];
funct.JmF = [];

funct.rF = [];
funct.rJ = [];
funct.rFmF = [];
funct.rJmF = [];

v = materials.v;
method = solParams.method;

if strcmp(method,'reducedNR') || strcmp(method,'hybridReducedNR')
    rootFindingMethod = 'reducedSystem';
else
    rootFindingMethod = 'fullSystem';
end


if (length(v)-1) == length(intF)
    n = length(v);  %Number of layers
    nElSub = size(source,1);
    Tx = target(1);
    Ty = target(2);
    
    syms('x',[(2*n)-1,1]);
    F = sym('F',[(2*n)-1,nElSub]);
    J = sym('J',[(2*n)-1,(2*n)-1,nElSub]);
    J(:,:,:) = sym(0);
    
    dIntF = sym('diffIntF',[(n-1) 1]);
    for i=1:(n-1)
        dIntF(i) = diff(intF(i),x(i));
    end
    
    for j=1:nElSub
        Sx = source(j,1);
        Sy = source(j,2);
        
        if n==1
            F(1,j) = ((Sy - Ty)*tan(x(n))) - Tx + Sx;
            J(1,n,j) = (Sy - Ty)*(sec(x(n))^2);
        else
            F(1,j) = ((Sy - (intF(1)))*tan(x(n))) - x(1) + Sx;
            J(1,1,j) = -1 -(tan(x(n))*dIntF(1));
            J(1,n,j) = (Sy - (intF(1)))*(sec(x(n))^2);
            
            for i=2:n-1
                F(i,j) = ((intF(i-1) - intF(i))*tan(x(n+i-1))) - x(i) + x(i-1);
                for k=1:((2*n)-1)
                    J(i,k,j) = diff(F(i,j),x(k));
                end
            end
            
            F(n,j) = ((intF(n-1) - Ty)*tan(x((2*n)-1))) - Tx + x(n-1);
            for k=1:((2*n)-1)
                J(n,k,j) = diff(F(n,j),x(k));
            end
            
            % Relationships between angles
            for i=2:n
                F(n+i-1,j) = asin((v(i)/v(i-1))*sin(x(n+i-2)-atan(dIntF(i-1))))+atan(dIntF(i-1)) - x(n+i-1);
                for k=1:((2*n)-1)
                    J(n+i-1,k,j) = diff(F(n+i-1,j),x(k));
                end
            end
        end
    end
    
    FmF = matlabFunction(F,'Vars',{x});
    JmF = matlabFunction(J,'Vars',{x});
    %             dIntFmF = matlabFunction(dIntF,'Vars',{x});
    
    funct.dIntF = dIntF;
    funct.F = F;
    funct.J = J;
    funct.FmF = FmF;
    funct.JmF = JmF;
    
    if strcmp(rootFindingMethod,'reducedSystem')
        syms('x',[n,1]);
        forwardAngles = sym('forwardAngles',[n,1]);
        forwardAngles(1) = x(n);
        F = sym('F',[n,nElSub]);
        J = sym('J',[n,n,nElSub]);
        J(:,:,:) = sym(0);
        
        dIntF = sym('diffIntF',[(n-1) 1]);
        for i=1:(n-1)
            dIntF(i) = diff(intF(i),x(i));
        end
        
        for j=1:nElSub
            Sx = source(j,1);
            Sy = source(j,2);
            
            if n==1
                F(1,j) = ((Sy - Ty)*tan(x(n))) - Tx + Sx;
                J(1,n,j) = (Sy - Ty)*(sec(x(n))^2);
            else
                F(1,j) = ((Sy - (intF(1)))*tan(x(n))) - x(1) + Sx;
                J(1,1,j) = -1 -(tan(x(n))*dIntF(1));
                J(1,n,j) = (Sy - (intF(1)))*(sec(x(n))^2);
                
                % Compute symbolic forward angles
                for i=2:n
                    forwardAngles(i) = asin((v(i)/v(i-1))*sin(forwardAngles(i-1)-atan(dIntF(i-1))))+atan(dIntF(i-1));
                end
                
                for i=2:n-1
                    F(i,j) = ((intF(i-1) - intF(i))*tan(forwardAngles(i))) - x(i) + x(i-1);
                    for k=1:n
                        J(i,k,j) = diff(F(i,j),x(k));
                    end
                end
                
                F(n,j) = ((intF(n-1) - Ty)*tan(forwardAngles(n))) - Tx + x(n-1);
                for k=1:n
                    J(n,k,j) = diff(F(n,j),x(k));
                end
            end
        end
        
        FmF = matlabFunction(F,'Vars',{x});
        JmF = matlabFunction(J,'Vars',{x});
        %             dIntFmF = matlabFunction(dIntF,'Vars',{x});
        
        funct.rF = F;
        funct.rJ = J;
        funct.rFmF = FmF;
        funct.rJmF = JmF;
    end
else
    F = nan;
    J = nan;
    funct.F = F;
    funct.J = J;
    if strcmp(rootFindingMethod,'reducedSystem')
        funct.rF = F;
        funct.rJ = J;
    end
    disp('Error! - Number of layer velocities is not equal to the number of interface functions plus 1.')
end

%------------- END CODE --------------

end

