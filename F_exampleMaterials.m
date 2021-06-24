function materials = F_exampleMaterials(n)
%materials = F_exampleMaterials(n) generates the material properties for n
%layers. This is just an ausiliary function to test the execution of the
%algorithms, relative to the paper titled: "Solving Ultrasonic Ray Tracing
%in Parts with Multiple Material Layers through Root-Finding Methods", by
%C. Mineo, D. Cerniglia, and E. Mohseni.
%
%   Inputs:
%       n - integer number giving the number of layers
%
%   Outputs:
%       materials - structured array
%           materials.v = [n x 1]  % wave propation speed in the material [mm/us]
%           materials.color = [n x 3]  % RGB color components for each material
%           materials.color = cell(n x 1)  % material names
%
% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------

materials = [];
materials.v = [1480,...  %Water 20°C
               2730,...  %Acrylic (Perspex)
               3320,...  %Tin
               4660,...  %Copper
               5850,...  %Steel 4340
               4430,...  %Brass
               6320];    %Aluminum

materials.color = [0.95 0.95 0.95;
                    0.90 0.90 0.80;
                    0.70 0.90 0.70;
                    0.90 0.80 0.80;
                    0.75 0.75 0.75;
                    0.80 0.70 0.60;
                    0.80 0.85 0.95];

materials.label = {'Water 20°C';
                    'Acrylic (Perspex)';
                    'Tin';
                    'Copper';
                    'Steel 4340';
                    'Brass';
                    'Aluminum'};

if n==1
    materials.v = materials.v(1);
    materials.color = materials.color(1,:);
    materials.label = materials.label{1};
elseif n==2
    materials.v = [materials.v(1);materials.v(end)];
    materials.color = [materials.color(1,:);materials.color(end,:)];
    materials.label = {materials.label{1};materials.label{end}};
elseif n==3
    materials.v = [materials.v(1);materials.v(3);materials.v(end)];
    materials.color = [materials.color(1,:);materials.color(3,:);materials.color(end,:)];
    materials.label = {materials.label{1};materials.label{3};materials.label{end}};
elseif n==4
    materials.v = [materials.v(1);materials.v(2);materials.v(4);materials.v(end)];
    materials.color = [materials.color(1,:);materials.color(2,:);materials.color(4,:);materials.color(end,:)];
    materials.label = {materials.label{1};materials.label{2};materials.label{4};materials.label{end}};
elseif n==5
    materials.v = [materials.v(1);materials.v(2);materials.v(4);materials.v(5);materials.v(end)];
    materials.color = [materials.color(1,:);materials.color(2,:);materials.color(4,:);materials.color(5,:);materials.color(end,:)];
    materials.label = {materials.label{1};materials.label{2};materials.label{4};materials.label{5};materials.label{end}};
elseif n==6
    materials.v = [materials.v(1);materials.v(2);materials.v(3);materials.v(4);materials.v(5);materials.v(end)];
    materials.color = [materials.color(1,:);materials.color(2,:);materials.color(3,:);materials.color(4,:);materials.color(5,:);materials.color(end,:)];
    materials.label = {materials.label{1};materials.label{2};materials.label{3};materials.label{4};materials.label{5};materials.label{end}};
elseif n==7
    materials.v = [materials.v(1);materials.v(2);materials.v(3);materials.v(4);materials.v(5);materials.v(6);materials.v(end)];
    materials.color = [materials.color(1,:);materials.color(2,:);materials.color(3,:);materials.color(4,:);materials.color(5,:);materials.color(6,:);materials.color(end,:)];
    materials.label = {materials.label{1};materials.label{2};materials.label{3};materials.label{4};materials.label{5};materials.label{6};materials.label{end}};
elseif n>7
    v1 = repmat(materials.v(2:end-1)',100,1);
    color1 = repmat(materials.color(2:end-1,:),100,1);
    ind = repmat(linspace(2,length(materials.v)-1,length(materials.v)-2)', 100, 1);
    label1 = materials.label(ind);

    v1 = v1(1:n-2);
    color1 = color1(1:n-2,:);
    ind = ind(1:n-2);
    
    materials.v = [materials.v(1);v1;materials.v(end)];
    materials.color = [materials.color(1,:);color1;materials.color(end,:)];
    materials.label = [materials.label(1);label1(ind);materials.label(end)];
end

%------------- END CODE --------------

end

