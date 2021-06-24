function [dIn,dOut,dTime,incPts,tn] = F_iterationDeviation(co,cn,fo,fn,to,probVars,intFmF,s,t,v)
%[dIn,dOut,dTime,incPts,tn] = F_iterationDeviation(co,cn,fo,fn,to,probVars,intFmF,s,t,v)
%computes the deviation between the last iteration of an iterative root-finding
%algorithms from the previous one
%
%   Inputs:
%       co - old test value
%       cn - new test value
%       fo - deviation relative to the old test value
%       fn - deviation relative to the new test value
%       to - propagation time of the ray at the old test value
%       probVars - problem variables
%       intFmF - symbolic matlab functions of the interface functions
%       s = [Sx,Sy] source coordinates
%       t = [Tx,Ty] target coordinates
%       v = propagation speed of all material layers
%
%   Outputs:
%       dIn - deviation between the input test values (dIn = cn-co)
%       dOut - deviation between the output test values (dOut = fn-fo)
%       dTime - deviation between propagation times
%       incPts - incidence points of the new ray path
%       tn - propagation time of the ray at the new test value
%
% Author: Carmelo Mineo
% Department of Engineering, University of Palermo, Viale delle Scienze,
% Edificio 8, 90128 Palermo, Italy.
% email: carmelo.mineo01@unipa.it
% Website: http://www.unipa.it
% June 2021; Last revision: 24-June-2021
% Tested with: Matlab 2020b


%------------- BEGIN CODE --------------

n = length(v);  %Number of layers

incPts = zeros(n+1,2);
incPts(1,:) = s;
for i=1:n-1
    incPts(i+1,:) = [probVars(i) intFmF{i}(probVars)];
end

if length(cn)>1
    cn = cn(n);
end

if length(co)>1
    co = co(n);
end

dIn = cn-co;
dOut = fn-fo;

incPts(n+1,:) = t;
dists = sqrt(sum(diff(incPts).^2,2));
travTimes = dists./v;
tn = sum(travTimes)*(10^6);
dTime = to - tn;

tx = t(1) + fn;
incPts(n+1,:) = [tx,t(2)];

%------------- END CODE --------------

end

