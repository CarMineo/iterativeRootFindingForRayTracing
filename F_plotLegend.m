function hl = F_plotLegend(hp,materials,fontsize,location)
%hl = F_plotLegend(hp,materials,fontsize,location) plots the legend of the
%display

% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------


ch = get(hp,'children');
hLegend = [];
legLabels = {};
hf = findobj(ch,'Type','line','Color',[0.4 0.4 0.4]);
if ~isempty(hf)
    hLegend(end+1) = hf(1);
    legLabels{end+1} = 'Initial guess';
end

hf = findobj(ch,'Type','line','Color',[0.3 0.3 1]);
if ~isempty(hf)
    hLegend(end+1) = hf(1);
    legLabels{end+1} = 'Bisection iteration';
end

hf = findobj(ch,'Type','line','Color',[1 0.2 0.2]);
if ~isempty(hf)
    hLegend(end+1) = hf(1);
    legLabels{end+1} = 'N-R iteration';
end

hf = findobj(ch,'Type','line','Color',[0.8 0.2 0.8]);
if ~isempty(hf)
    hLegend(end+1) = hf(1);
    legLabels{end+1} = 'Reduced N-R iteration';
end

hf = findobj(ch,'Type','line','Color',[0 0 0]);
if ~isempty(hf)
    hLegend(end+1) = hf(1);
    legLabels{end+1} = 'Final result';
end

[uniqueColors,ui] = unique(materials.color,'rows','stable');
uniqueLabels = materials.label(ui);

nu = length(ui);
for i=1:nu
    hf = findobj(ch,'Type','patch','faceColor',uniqueColors(i,:));
    if ~isempty(hf)
        hLegend(end+1) = hf(1);
        legLabels{end+1} = uniqueLabels{i};
    end
end

hl = legend(hLegend,legLabels,'FontSize',fontsize,'Location',location);

%------------- END CODE --------------
