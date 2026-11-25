function C = capacitance(G, epsr, mask, Vfix, opts)
%CAPACITANCE Two-terminal capacitance per unit depth by stored energy.
%   C=EM.DERIVE.CAPACITANCE(G,EPSR,MASK,VFIX,OPTS) returns F/m from
%   C=2*energyE(G,V,epsr)/deltaV^2 after a zero-charge Poisson solve.
%   MASK/VFIX must specify exactly two distinct conductor potentials.
%   A common offset and arbitrary nonzero voltage difference are allowed.
%   OPTS is passed to poisson2d; harmonic faces are required so the energy
%   and solver represent the same material model.
%
%   Example:
%       G=em.solve.grid2d([0 1],[0 1],21,21);
%       m=false(21); m([1 end],:)=true; v=zeros(21); v(end,:)=1;
%       C=em.derive.capacitance(G,1,m,v,struct('method','direct'));
%
%   See also EM.DERIVE.ENERGYE, EM.DERIVE.RESISTANCE.
if nargin<5||isempty(opts), opts=struct(); end
if ~isstruct(opts)||~isscalar(opts), error('em:derive:capacitance:Opts','opts must be a scalar struct.'); end
if isfield(opts,'faceMean')&&~strcmpi(opts.faceMean,'harmonic'), error('em:derive:capacitance:Mean','opts.faceMean must be harmonic for energy consistency.'); end
if ~islogical(mask)||~isnumeric(Vfix)||~isequal(size(mask),size(Vfix))||~any(mask(:)), error('em:derive:capacitance:Mask','mask and Vfix must be matching NyxNx arrays.'); end
levels=unique(Vfix(mask));
if numel(levels)~=2||~isreal(levels)||any(~isfinite(levels)), error('em:derive:capacitance:Voltage','Vfix(mask) must contain two finite distinct conductor levels.'); end
V=em.solve.poisson2d(G,0,epsr,mask,Vfix,opts);
C=2*em.derive.energyE(G,V,epsr)/diff(levels)^2;
end
