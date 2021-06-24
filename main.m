close all; clear; clc;

%This matlab script demonstrates the execution of the algorithms, relative
%to the paper titled: "Solving Ultrasonic Ray Tracing in Parts with
%Multiple Material Layers through Root-Finding Methods", by C. Mineo,
%D. Cerniglia, and E. Mohseni.
%
% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         INPUT PARAMETERS                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = [0,0];       % Wave source point coordinates [x[mm], y[mm]]
t = [20,-30];    % Wave target point coordinates [x[mm], y[mm]]
n = 4;           % Number of material layers

% Phased-Array parameters (set nElSub=1 for a single element source)
PAangle = 0;    % Angle the PA probe forms with the horizontal direction [degrees]
nElSub = 1;     % Number of elements in the subaperture
PApitch = 1;    % Pitch of the Phased Array probe [mm]
elWidth = 0.5;  % Width of each phased array element [mm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of materials (velocity [mm/us], plot color and label) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: In this example a set of velocities (for seven different materials)
%       is defined and the following code selects how to arrange the order
%       of the layers, to fill "n" number of layers. You can replace these
%       velocities with other values, omit the use of "F_exampleMaterials"
%       and define your preferred materials for the sequence of n layers.

% materials.v = [n x 1]  % wave propation speed in the material [mm/us]
% materials.color = [n x 3]  % RGB color components for each material
% materials.color = cell(n x 1)  % material names
materials = F_exampleMaterials(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS OF INTERFACES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: The functions describing the interfaces can be flat (constant),
% polynomial (of any degree), trigonometric, or any other type of
% analytical functions. Important requirement is that the functions are
% derivable and the derivative are continuous.

syms('x',[n,1]);                % Initialise a vector of symbolic variables
intF = sym('intF',[(n-1) 1]);   % Initialise a vector of symbolic functions
intFmF = cell(n-1,1);           % Initialise a cell array to store matlab functions

% In this example trigonometric functions are defined by the following
% code. You can replace this code with custom code for the definition of
% your preferred interface functions.
if n>1
    for i=1:n-1
        amp = ((s(2)-t(2))/(4*n));
        intF(i) = amp*sind(x(i)*10 + ((i-1)*45)) - ((s(2)-t(2))/n)*i;
        intFmF{i} = matlabFunction(intF(i),'Vars',{x});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR ITERATIVE SOLVING METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solParams = [];

%solParams.method = 'bisection';  
%solParams.method = 'bisectionWithNestedNR';
%solParams.method = 'NR';
%solParams.method = 'reducedNR';
solParams.method = 'hybridNR';
%solParams.method = 'hybridReducedNR';

solParams.numTol = eps('single');     % Maximum numerical tolerance
solParams.stopMode = 'dt';            % Stopping condition
solParams.stopTol = 0.001; %[us] (Tolerance set to 1ns for this example)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLAGS FOR DISPLAY AND SAVING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iterPlotting = true;
showLegend = true;
saveFigure = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     COMPUTATION OF SOLUTION                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% INITIALISE FIGURE %
%%%%%%%%%%%%%%%%%%%%%

hFig=figure(1); set(hFig, 'Position', get(0, 'Screensize'));
[~,~,xmin,xmax] = F_plotLayers(intF,materials,s,t,PAangle,PApitch,elWidth,nElSub);
set(gca,'fontsize',14); xlabel('x [mm]'); ylabel('y [mm]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SYMBOLIC FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

funct = F_genFunctions(intF,materials,s,t,solParams);

%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE INITIAL GUESS %
%%%%%%%%%%%%%%%%%%%%%%%%%

initGuess = F_findInitialGuess(funct,materials,s,t,xmin,xmax,solParams);
[~,incPts,solPerform] = F_solveRayTracing(funct,initGuess,materials,s,t,xmin,xmax,solParams,iterPlotting);

%%%%%%%%%%%%%%%%%%%%%
% PLOT FINAL RESULT %
%%%%%%%%%%%%%%%%%%%%%
plot(incPts(:,1),incPts(:,2),'.-k','markersize',25,'linewidth',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT LEGEND (IF REQUESTED) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if showLegend
    F_plotLegend(gca,materials,16,'northeastoutside');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE FIGURE (IF REQUESTED) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveFigure
    savefig(hFig,[solParams.method '_' num2str(n) '.fig']);
    print([solParams.method '_' num2str(n) '.tiff'],'-dtiff','-r300');
end