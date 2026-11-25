function W = energyE(G, V, epsr)
%ENERGYE Electric energy per unit depth from the discrete grid gradient.
%   W=EM.DERIVE.ENERGYE(G,V,EPSR) returns J/m, using one half eps0
%   times the sum of harmonic face epsr times squared edge voltage drops.
%   V is NyxNx volts; EPSR is positive scalar or NyxNx. Face gradients
%   and dual cell areas give trapezoidal half weights at outer boundaries.
%   This is the same conservative harmonic stencil used by poisson2d.
%   Interfaces are halfway between nodes; conductor interior edges have
%   zero gradient. No interpolant or pointwise central gradient is used.
%
%   Example:
%       G=em.solve.grid2d([0 1],[0 1],21,21);
%       W=em.derive.energyE(G,G.X,1); % eps0/2 J/m
%
%   See also EM.DERIVE.CAPACITANCE, EM.SOLVE.POISSON2D.
if ~isstruct(G)||~all(isfield(G,{'X','Y','hx','hy'}))||~isscalar(G.hx)||~isscalar(G.hy)||~isfinite(G.hx)||~isfinite(G.hy)||G.hx<=0||G.hy<=0
    error('em:derive:energyE:Grid','G must come from grid2d.');
end
sz=size(G.X);
if numel(sz)~=2||min(sz)<2||~isequal(size(G.Y),sz), error('em:derive:energyE:Grid','G must contain NyxNx arrays.'); end
if ~isnumeric(V)||~isreal(V)||~isequal(size(V),sz)||any(~isfinite(V(:))), error('em:derive:energyE:V','V must be finite real NyxNx.'); end
if isnumeric(epsr)&&isscalar(epsr), epsr=repmat(epsr,sz); end
if ~isnumeric(epsr)||~isreal(epsr)||~isequal(size(epsr),sz)||any(~isfinite(epsr(:)))||any(epsr(:)<=0), error('em:derive:energyE:Epsr','epsr must be positive finite scalar or NyxNx.'); end
wx=ones(1,sz(2)); wx([1 end])=0.5; wy=ones(sz(1),1); wy([1 end])=0.5;
ex=2./(1./epsr(:,1:end-1)+1./epsr(:,2:end));
ey=2./(1./epsr(1:end-1,:)+1./epsr(2:end,:));
a=ex.*diff(V,1,2).^2.*wy*(G.hy/G.hx);
b=ey.*diff(V,1,1).^2.*wx*(G.hx/G.hy);
W=0.5*em.const.eps0()*(sum(a(:))+sum(b(:)));
end
