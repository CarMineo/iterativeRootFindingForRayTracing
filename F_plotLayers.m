function [handle,centreOfElements,xmin,xmax] = F_plotLayers(intF,materials,s,t,PAangle,PApitch,elWidth,nElSub)
%[handle,centreOfElements,xmin,xmax] = F_plotLayers(intF,materials,s,t,PAangle,PApitch,elWidth,nElSub)

% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------


v = materials.v;

%%%% PLOT LAYERS %%%
lPA = (nElSub - 1)*PApitch;     % Length of subaperture in Phased Array probe
n = size(v,1);
syms('x',[n,1]);

subapertureCentre = 0;
offset = 0;

%%% PLOT SYMBOLIC REPRESENTATION OF PHASED ARRAY SUBAPERTURE %%%
startOfElements = (-lPA/2-elWidth/2):PApitch:(lPA/2-elWidth/2);
startOfElements = startOfElements + offset;

BLX = startOfElements*cosd(PAangle)+subapertureCentre;
BLY = startOfElements*sind(PAangle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE FOCUSING PATHS FOR ALL ELEMENTS OF THE SUBAPERTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centreOfElementsX = BLX + (elWidth/2)*cosd(PAangle);
centreOfElementsY = BLY + (elWidth/2)*sind(PAangle);
centreOfElements = zeros(nElSub,2);
centreOfElements(:,1) = centreOfElementsX;
centreOfElements(:,2) = centreOfElementsY;

if nElSub==1
    centreOfElements(:,1) = centreOfElements(:,1)-centreOfElements(1,1) + s(1);
    centreOfElements(:,2) = centreOfElements(:,2)-centreOfElements(1,2) + s(2);
else
    centreOfElements(:,1) = centreOfElements(:,1) + s(1);
    centreOfElements(:,2) = centreOfElements(:,2) + s(2);
end

plot(s(1),s(2),'.k'); hold on;

xmin = min(min(centreOfElements(:,1)),t(1));
xmax = max(max(centreOfElements(:,1)),t(1));
xmin = xmin - 0.1*(xmax - xmin);
xmax = xmax + 0.1*(xmax - xmin);

ymin = min(min(centreOfElements(:,2)),t(2));
ymax = max(max(centreOfElements(:,2)),t(2));
ymin = ymin - 0.11*(ymax - ymin);
ymax = ymax + 0.1*(ymax - ymin);




xArray = (xmin-1):0.5:(xmax+1);
for i=1:n
    if n==1
        h = ymax + 0.1*(ymax - ymin);
        yt = h*ones(size(xArray));
        yb = ones(size(xArray)).*(ymin - 0.1*(ymax - ymin));
        xp = [xArray';fliplr(xArray)';xArray(1)];
        yp = [yb';yt';yb(1)];
    else
        if i==1
            yb = double(subs(intF(i),x(i),xArray));
            if length(yb)==1
                yb = ones(size(xArray)).*yb;
            end
            h = ymax + 0.1*(ymax - ymin);
            yt = h*ones(size(xArray));
            xp = [xArray';xArray(end);xArray(1);xArray(1)];
            yp = [yb';yt(end);yt(1);yb(1)];
        elseif i==n
            yt = double(subs(intF(i-1),x(i-1),fliplr(xArray)));
            if length(yt)==1
                yt = ones(size(xArray)).*yt;
            end
            yb = ones(size(xArray)).*(ymin - 0.1*(ymax - ymin));
            xp = [xArray';fliplr(xArray)';xArray(1)];
            yp = [yb';yt';yb(1)];
        else
            yt = double(subs(intF(i-1),x(i-1),fliplr(xArray)));
            if length(yt)==1
                yt = ones(size(xArray)).*yt;
            end
            yb = double(subs(intF(i),x(i),xArray));
            if length(yb)==1
                yb = ones(size(xArray)).*yb;
            end
            xp = [xArray';fliplr(xArray)';xArray(1)];
            yp = [yb';yt';yb(1)];
        end
    end
    col = materials.color(i,:);
    fill(xp,yp,col,'linewidth',1.5);hold on;
end

%%% INDICATE SOURCE POINT WITH A BLACK DOT %%%
for i=1:nElSub
    plot(centreOfElements(i,1), centreOfElements(i,2),'.k','markersize',27); hold on;
    if nElSub==1
        text(centreOfElements(i,1)-0.7, centreOfElements(i,2)+1.6,'S','FontSize',32,'FontWeight','bold');
    else
        if i==1
            plot(s(1), s(2),'*r','markersize',15,'linewidth',2);
            text(s(1)-0.3, s(2)+1,'O','FontSize',15,'FontWeight','bold');
        end
    end
end

%%% INDICATE TOP OF FIRST LAYER WITH DOTTED LINE %%%
plot([(xmin-1) (xmax+1)],[s(2) s(2)],'--k','linewidth',1.5); hold on;

%%% INDICATE BOTTOM OF LAST LAYER WITH DOTTED LINE %%%
plot([(xmin-1) (xmax+1)],[t(2) t(2)],'--k','linewidth',1.5);

%%% INDICATE TARGET POINT WITH A BLACK DOT %%%
plot(t(1), t(2),'.k','markersize',27);
text(t(1)-0.3, t(2)-1.4,'T','FontSize',32,'FontWeight','bold');

%pbaspect([1 1 1]);
axis equal;

%%% LIMIT WINDOW IN X-DIRECTION %%%
xlim([xmin xmax]);

%%% LIMIT WINDOW IN Y-DIRECTION %%%
ylim([ymin ymax]);

handle = gca;

%------------- END CODE --------------

