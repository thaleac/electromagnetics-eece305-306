function R = resistance(G, sigma, mask, Vfix, opts)
%RESISTANCE Two-terminal resistance times unit depth, in ohm m.
%   R=EM.DERIVE.RESISTANCE(G,SIGMA,MASK,VFIX,OPTS) solves div(sigma
%   grad V)=0. SIGMA is positive scalar or NyxNx S/m. Current per depth
%   leaving the high-potential conductor is summed over every incident
%   grid face using -sigma*dV/dn times face length. R=deltaV/Iprime.
%   MASK/VFIX specify exactly two distinct terminal potentials. Unmasked
%   outer faces are insulating. OPTS is passed to poisson2d, including
%   faceMean (default harmonic). A depth d gives actual resistance R/d.
%
%   Example:
%       G=em.solve.grid2d([0 1],[0 1],21,21);
%       m=false(21); m(:,[1 end])=true; v=zeros(21); v(:,1)=1;
%       R=em.derive.resistance(G,100,m,v,struct('method','direct'));
%
%   See also EM.DERIVE.CAPACITANCE, EM.SOLVE.POISSON2D.
if nargin<5||isempty(opts), opts=struct(); end
if ~islogical(mask)||~isnumeric(Vfix)||~isequal(size(mask),size(Vfix))||~any(mask(:)), error('em:derive:resistance:Mask','mask and Vfix must be matching NyxNx arrays.'); end
levels=unique(Vfix(mask));
if numel(levels)~=2||~isreal(levels)||any(~isfinite(levels)), error('em:derive:resistance:Voltage','Vfix(mask) must contain two finite distinct terminal levels.'); end
V=em.solve.poisson2d(G,0,sigma,mask,Vfix,opts);
if isscalar(sigma), sigma=repmat(sigma,size(V)); end
high=mask&Vfix==levels(2); wx=ones(1,size(V,2));wx([1 end])=0.5;
wy=ones(size(V,1),1);wy([1 end])=0.5;
kind='harmonic';if isfield(opts,'faceMean'),kind=opts.faceMean;end
sx=face(sigma(:,1:end-1),sigma(:,2:end),kind);
sy=face(sigma(1:end-1,:),sigma(2:end,:),kind);
ix=sx.*(-diff(V,1,2)).*(-diff(double(high),1,2)).*wy*(G.hy/G.hx);
iy=sy.*(-diff(V,1,1)).*(-diff(double(high),1,1)).*wx*(G.hx/G.hy);
I=sum(ix(:))+sum(iy(:));
if ~(isfinite(I)&&I>0), error('em:derive:resistance:Current','sigma and terminal geometry must give positive finite current.'); end
R=diff(levels)/I;
end
function f=face(a,b,kind)
if strcmpi(kind,'harmonic'),f=2./(1./a+1./b);else,f=(a+b)/2;end
end
