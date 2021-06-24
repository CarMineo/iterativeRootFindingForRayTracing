function [probVars,incPts,solPerform] = F_solveRayTracing(funct,initGuess,materials,s,t,xmin,xmax,solParams,interPlotting)
%[probVars,incPts,solPerform] = F_solveRayTracing(funct,initGuess,materials,s,t,xmin,xmax,solParams,interPlotting)
%solves the ray tracing problem through one of the iterative root-finding
%algorithms, investigated in the paper titled: "Solving Ultrasonic Ray
%Tracing in Parts with Multiple Material Layers through Root-Finding
%Methods", by C. Mineo, D. Cerniglia, and E. Mohseni.
%
%   Inputs:
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
%       initGuess - structured array
%           initGuess.a = a real number giving the left extremity of the
%                         search interval for the bisection-based methods;
%           initGuess.b = a real number giving the right extremity of the
%                         search interval for the bisection-based methods;
%           initGuess.c = a real number giving the central value of the
%                         search interval for the bisection-based methods;
%           initGuess.fa = horizontal deviation in initGuess.a;
%           initGuess.fc = horizontal deviation in initGuess.c;
%       materials - structured array
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
%       interPlotting - a boolean flag to specify if interim ray paths are
%                       plotted (interPlotting = true) or not (interPlotting = false)
%
%   Outputs:
%       probVar - problem variables
%       incPts - coordinates of incidence points between the ray path and
%                the interfaces 
%       solPerform - structured array containing the solution performance
%           solPerform.dx - deviation between input variables
%           solPerform.dy - difference between output deviations
%           solPerform.dt - difference between ray travelling time
%           solPerform.solTime - total elapsed time to reach the solution;
%           solPerform.bisIterTime - average computation time of bisection iterations;
%           solPerform.NRiterTime - average computation time of N-R iterations;
%           solPerform.nBisection - number of bisection iterations;
%           solPerform.nNewton - number of N-R iterations;
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

intF = funct.intF;

F = funct.F;
J = funct.J;
FmF = funct.FmF;
JmF = funct.JmF;

rFmF = funct.rFmF;
rJmF = funct.rJmF;

a = initGuess.a;
b = initGuess.b;
c = initGuess.c;
fa = initGuess.fa;
fc = initGuess.fc;

numTol = solParams.numTol;
stopTol = solParams.stopTol;
stopMode = solParams.stopMode;
method = solParams.method;

n = length(v);  %Number of layers

syms('x',[(2*n)-1,1]);
intFmF = cell(n-1,1);
for i=1:n-1
    intFmF{i} = matlabFunction(intF(i),'Vars',{x});
end

nBisection = 0;
nNewton = 0;
bisIterTime = 0;
NRiterTime = 0;

switch method
    case 'bisection'
        % PURE BISECTION METHOD
        newtonMode = false;
        
        co = inf;
        fo = inf;
        cn = c(n);
        to = inf;
        probVars = c;
        
        [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc,to,probVars,intFmF,s,t,v);
        switch stopMode
            case 'dx'
                dev = abs(dx);
            case 'dy'
                dev = abs(dy);
            case 'dt'
                dev = abs(dt);
        end
        
        if interPlotting
            plot(incPts(:,1),incPts(:,2),'.:','color',[0.4 0.4 0.4],'markersize',20,'linewidth',5);
            if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                text(incPts(end,1)-0.3, incPts(end,2)-1.6,'T_0','color',[0.2 0.2 0.2],'FontSize',24,'FontWeight','bold');
            end
            drawnow;
        end
        
        solStart = tic;
        
        iterationsCounter = 0;
        while dev > stopTol
            iterStart = tic;
            
            iterationsCounter = iterationsCounter + 1;
            nBisection = nBisection + 1;
            
            if sign(fc)~=sign(fa)
                b = cn;
            else
                a = cn;
                fa = fc;
            end
            
            fo = fc;
            to = tn;
            co = cn;
            cn = (a+b)/2;
            [probVars,fc] = F_propagateRay(intF,F,J,s,cn,xmin,xmax,numTol,newtonMode);
            [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc,to,probVars,intFmF,s,t,v);
            
            bisIterTime = bisIterTime + toc(iterStart);
            
            switch stopMode
                case 'dx'
                    dev = abs(dx);
                case 'dy'
                    dev = abs(dy);
                case 'dt'
                    dev = abs(dt);
            end
            
            if interPlotting
                plot(incPts(:,1),incPts(:,2),'.:','color',[0.3 0.3 1],'markersize',20,'linewidth',5);
                if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                    text(incPts(end,1)-0.3, incPts(end,2)-1.6,['T_' num2str(iterationsCounter)],'color',[0.1 0.1 1],'FontSize',24,'FontWeight','bold');
                end
                drawnow;
            end
        end
        solTime = toc(solStart);
        bisIterTime = bisIterTime/nBisection;
        
    case 'bisectionWithNestedNR'
        % HYBRID BISECTION METHOD
        newtonMode = true;
        
        co = inf;
        fo = inf;
        cn = c(n);
        to = inf;
        probVars = c;
        
        [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc,to,probVars,intFmF,s,t,v);
        switch stopMode
            case 'dx'
                dev = abs(dx);
            case 'dy'
                dev = abs(dy);
            case 'dt'
                dev = abs(dt);
        end
        
        if interPlotting
            plot(incPts(:,1),incPts(:,2),'.:','color',[0.4 0.4 0.4],'markersize',20,'linewidth',5);
            if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                text(incPts(end,1)-0.3, incPts(end,2)-1.6,'T_0','color',[0.2 0.2 0.2],'FontSize',24,'FontWeight','bold');
            end
            drawnow;
        end
        
        solStart = tic;
        
        iterationsCounter = 0;
        while dev > stopTol
            iterStart = tic;
            
            iterationsCounter = iterationsCounter + 1;
            nBisection = nBisection + 1;
            
            if sign(fc)~=sign(fa)
                b = cn;
            else
                a = cn;
                fa = fc;
            end
            
            fo = fc;
            to = tn;
            co = cn;
            cn = (a+b)/2;
            [probVars,fc] = F_propagateRay(intF,F,J,s,cn,xmin,xmax,numTol,newtonMode);
            [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc,to,probVars,intFmF,s,t,v);
            
            bisIterTime = bisIterTime + toc(iterStart);
            
            switch stopMode
                case 'dx'
                    dev = abs(dx);
                case 'dy'
                    dev = abs(dy);
                case 'dt'
                    dev = abs(dt);
            end
            
            if interPlotting
                plot(incPts(:,1),incPts(:,2),'.:','color',[0.3 0.3 1],'markersize',20,'linewidth',5);
                if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                    text(incPts(end,1)-0.3, incPts(end,2)-1.6,['T_' num2str(iterationsCounter)],'color',[0.1 0.1 1],'FontSize',24,'FontWeight','bold');
                end
                drawnow;
            end
        end
        solTime = toc(solStart);
        bisIterTime = bisIterTime/iterationsCounter;

    case 'NR'
        % FULL NEWTON-RAPHSON METHOD
        newtonMode = false;
        
        co = inf;
        fo = inf;
        cn = c;
        fc = FmF(cn);
        to = inf;
        probVars = c;
        
        [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
        switch stopMode
            case 'dx'
                dev = abs(dx);
            case 'dy'
                dev = abs(dy);
            case 'dt'
                dev = abs(dt);
        end
        
        if interPlotting
            plot(incPts(:,1),incPts(:,2),'.:','color',[0.4 0.4 0.4],'markersize',20,'linewidth',5);
            if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                text(incPts(end,1)-0.3, incPts(end,2)-1.6,'T_0','color',[0.2 0.2 0.2],'FontSize',24,'FontWeight','bold');
            end
            drawnow;
        end
        
        solStart = tic;
        
        iterationsCounter = 0;
        while dev > stopTol
            iterStart = tic;
            
            iterationsCounter = iterationsCounter + 1;
            nNewton = nNewton + 1;
            
            co = cn;
            fo = fc(n);
            to = tn;
            
            J = JmF(cn);
            cn = cn - (J\fc);
            fc = FmF(cn);
            
            ignore = false;
            if (cn(n)<a) || (cn(n)>b) || (~isreal(J)) || (~isreal(cn)) || (~isreal(fc))
                nNewton = nNewton - 1;
                ignore = true;
            end
            
            probVars = cn;
            
            [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
            
            if ignore==false
                NRiterTime = NRiterTime + toc(iterStart);
            end
            
            switch stopMode
                case 'dx'
                    dev = abs(dx);
                case 'dy'
                    dev = abs(dy);
                case 'dt'
                    dev = abs(dt);
            end
            
            if interPlotting
                plot(incPts(:,1),incPts(:,2),'.:','color',[1 0.2 0.2],'markersize',20,'linewidth',5);
                if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                    text(incPts(end,1)-0.3, incPts(end,2)-1.6,['T_' num2str(iterationsCounter)],'color',[1 0.1 0.1],'FontSize',24,'FontWeight','bold');
                end
                drawnow;
            end
        end
        
        solTime = toc(solStart);
        NRiterTime = NRiterTime/nNewton;
        
case 'hybridNR'
        % FULL NEWTON-RAPHSON METHOD
        newtonMode = false;
        
        co = inf;
        fo = inf;
        cn = c;
        fc = FmF(cn);
        to = inf;
        probVars = c;
        
        [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
        switch stopMode
            case 'dx'
                dev = abs(dx);
            case 'dy'
                dev = abs(dy);
            case 'dt'
                dev = abs(dt);
        end
        
        if interPlotting
            plot(incPts(:,1),incPts(:,2),'.:','color',[0.4 0.4 0.4],'markersize',20,'linewidth',5);
            if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                text(incPts(end,1)-0.3, incPts(end,2)-1.6,'T_0','color',[0.2 0.2 0.2],'FontSize',24,'FontWeight','bold');
            end
            drawnow;
        end
        
        solStart = tic;
        
        iterationsCounter = 0;
        while dev > stopTol
            iterStart = tic;
            
            iterationsCounter = iterationsCounter + 1;
            
            co = cn;
            fo = fc(n);
            to = tn;
            
            J = JmF(cn);
            cn = cn - (J\fc);
            fc = FmF(cn);
            
            if (cn(n)<a) || (cn(n)>b) || (~isreal(J)) || (~isreal(cn)) || (~isreal(fc))
                nBisection = nBisection + 1;
                if sign(fo)~=sign(fa)
                    b = co(n);
                else
                    a = co(n);
                    fa = fo;
                end
                
                cn = (a+b)/2;
                [cn,~] = F_propagateRay(intF,F,J,s,cn,xmin,xmax,numTol,newtonMode);
                
                fc = FmF(cn);
                probVars = cn;
                [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
                
                bisIterTime = bisIterTime + toc(iterStart);
                
                if interPlotting
                    plot(incPts(:,1),incPts(:,2),'.:','color',[0.3 0.3 1],'markersize',20,'linewidth',5);
                    if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                        text(incPts(end,1)-0.3, incPts(end,2)-1.6,['T_' num2str(iterationsCounter)],'color',[0.1 0.1 1],'FontSize',24,'FontWeight','bold');
                    end
                    drawnow;
                end
            else
                nNewton = nNewton + 1;
                probVars = cn;
                [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
                
                NRiterTime = NRiterTime + toc(iterStart);
                
                if interPlotting
                    plot(incPts(:,1),incPts(:,2),'.:','color',[1 0.2 0.2],'markersize',20,'linewidth',5);
                    if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                        text(incPts(end,1)-0.3, incPts(end,2)-1.6,['T_' num2str(iterationsCounter)],'color',[1 0.1 0.1],'FontSize',24,'FontWeight','bold');
                    end
                    drawnow;
                end
            end
            
            switch stopMode
                case 'dx'
                    dev = abs(dx);
                case 'dy'
                    dev = abs(dy);
                case 'dt'
                    dev = abs(dt);
            end
            
        end
        
        solTime = toc(solStart);
        bisIterTime = bisIterTime/nBisection;
        NRiterTime = NRiterTime/nNewton;
        
    case 'reducedNR'
        % FULL NEWTON-RAPHSON METHOD
        newtonMode = false;
        
        co = inf;
        fo = inf;
        cn = c(1:n);
        fc = rFmF(cn);
        to = inf;
        
        probVars = c;
        
        [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
        switch stopMode
            case 'dx'
                dev = abs(dx);
            case 'dy'
                dev = abs(dy);
            case 'dt'
                dev = abs(dt);
        end
        
        if interPlotting
            plot(incPts(:,1),incPts(:,2),'.:','color',[0.4 0.4 0.4],'markersize',20,'linewidth',5);
            if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                text(incPts(end,1)-0.3, incPts(end,2)-1.6,'T_0','color',[0.2 0.2 0.2],'FontSize',24,'FontWeight','bold');
            end
            drawnow;
        end
        
        solStart = tic;
        
        iterationsCounter = 0;
        while dev > stopTol
            iterStart = tic;
            
            iterationsCounter = iterationsCounter + 1;
            nNewton = nNewton + 1;
            
            co = cn;
            fo = fc(n);
            to = tn;
            
            rJ = rJmF(cn);
            cn = cn - (rJ\fc);
            fc = rFmF(cn);
            
            ignore = false;
            if (cn(n)<a) || (cn(n)>b) || (~isreal(rJ)) || (~isreal(cn)) || (~isreal(fc))
                nNewton = nNewton - 1;
                ignore = true;
            end
            
            probVars = cn;
            
            [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
            
            if ignore==false
                NRiterTime = NRiterTime + toc(iterStart);
            end
            
            switch stopMode
                case 'dx'
                    dev = abs(dx);
                case 'dy'
                    dev = abs(dy);
                case 'dt'
                    dev = abs(dt);
            end
            
            if interPlotting
                plot(incPts(:,1),incPts(:,2),'.:','color',[0.8 0.2 0.8],'markersize',20,'linewidth',5);
                if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                    text(incPts(end,1)-0.3, incPts(end,2)-1.6,['T_' num2str(iterationsCounter)],'color',[0.7 0.1 0.7],'FontSize',24,'FontWeight','bold');
                end
                drawnow;
            end
        end
        
        solTime = toc(solStart);
        NRiterTime = NRiterTime/nNewton;
        
case 'hybridReducedNR'
        % FULL NEWTON-RAPHSON METHOD
        newtonMode = false;
        
        co = inf;
        fo = inf;
        cn = c(1:n);
        fc = rFmF(cn);
        to = inf;
        probVars = c;
        
        [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
        
        switch stopMode
            case 'dx'
                dev = abs(dx);
            case 'dy'
                dev = abs(dy);
            case 'dt'
                dev = abs(dt);
        end
        
        if interPlotting
            plot(incPts(:,1),incPts(:,2),'.:','color',[0.4 0.4 0.4],'markersize',20,'linewidth',5);
            if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                text(incPts(end,1)-0.3, incPts(end,2)-1.6,'T_0','color',[0.2 0.2 0.2],'FontSize',24,'FontWeight','bold');
            end
            drawnow;
        end
        
        solStart = tic;
        
        iterationsCounter = 0;
        while dev > stopTol
            iterStart = tic;
            
            iterationsCounter = iterationsCounter + 1;
            
            co = cn;
            fo = fc(n);
            to = tn;
            
            rJ = rJmF(cn);
            cn = cn - (rJ\fc);
            fc = rFmF(cn);
            
            if (cn(n)<a) || (cn(n)>b) || (~isreal(rJ)) || (~isreal(cn)) || (~isreal(fc))
                nBisection = nBisection + 1;
                if sign(fo)~=sign(fa)
                    b = co(n);
                else
                    a = co(n);
                    fa = fo;
                end
                
                cn = (a+b)/2;
                [cn,~] = F_propagateRay(intF,F,J,s,cn,xmin,xmax,numTol,newtonMode);
                
                fc = rFmF(cn(1:n));
                probVars = cn;
                [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
                cn = cn(1:n);
                
                bisIterTime = bisIterTime + toc(iterStart);
                
                if interPlotting
                    plot(incPts(:,1),incPts(:,2),'.:','color',[0.3 0.3 1],'markersize',20,'linewidth',5);
                    if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                        text(incPts(end,1)-0.3, incPts(end,2)-1.6,['T_' num2str(iterationsCounter)],'color',[0.1 0.1 1],'FontSize',24,'FontWeight','bold');
                    end
                    drawnow;
                end
            else
                nNewton = nNewton + 1;
                probVars = cn;
                [dx,dy,dt,incPts,tn] = F_iterationDeviation(co,cn,fo,fc(n),to,probVars,intFmF,s,t,v);
                
                NRiterTime = NRiterTime + toc(iterStart);
                
                if interPlotting
                    plot(incPts(:,1),incPts(:,2),'.:','color',[0.8 0.2 0.8],'markersize',20,'linewidth',5);
                    if (abs(incPts(end,1)-t(1))>(0.1*(xmax-xmin))) && (incPts(end,1)>(xmin+(0.02*(xmax-xmin)))) && (incPts(end,1)<(xmax-(0.1*(xmax-xmin))))
                        text(incPts(end,1)-0.3, incPts(end,2)-1.6,['T_' num2str(iterationsCounter)],'color',[0.7 0.1 0.7],'FontSize',24,'FontWeight','bold');
                    end
                    drawnow;
                end
            end
            
            switch stopMode
                case 'dx'
                    dev = abs(dx);
                case 'dy'
                    dev = abs(dy);
                case 'dt'
                    dev = abs(dt);
            end
            
        end
        
        solTime = toc(solStart);
        bisIterTime = bisIterTime/nBisection;
        NRiterTime = NRiterTime/nNewton;
end

solPerform = [];
solPerform.dx = dx;
solPerform.dy = dy;
solPerform.dt = dt;
solPerform.solTime = solTime;
solPerform.bisIterTime = bisIterTime;
solPerform.NRiterTime = NRiterTime;
solPerform.nBisection = nBisection;
solPerform.nNewton = nNewton;

%------------- END CODE --------------