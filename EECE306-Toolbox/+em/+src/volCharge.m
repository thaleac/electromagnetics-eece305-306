function s = volCharge(rhoV, volume, uspan, vspan, wspan, Nu, Nv, Nw, rule)
%VOLCHARGE Construct a discretized volume-charge source.
%   S = EM.SRC.VOLCHARGE(RHOV, VOLUME, USPAN, VSPAN, WSPAN, NU, NV, NW,
%   RULE) accepts a scalar density in C/m^3 or a function handle RHOV(POS).
%   The default quadrature rule is 'midpoint'.
%
%   Example:
%       v = @(u,v,w) [u v w];
%       s = em.src.volCharge(1e-9,v,[0 1],[0 1],[0 1],10,10,10);
%
%   See also EM.QUAD.VOL, EM.SRC.LINECHARGE, EM.SRC.SURFCHARGE.

if nargin < 9 || isempty(rule), rule = 'midpoint'; end
[pos,w] = em.quad.vol(volume, uspan, vspan, wspan, Nu, Nv, Nw, rule);
q = evaluateDensity(rhoV, pos, 'rhoV');
s = struct('type','charge', 'pos',pos, 'w',w, 'q',q, ...
    'label','volume charge');
end

function q = evaluateDensity(rho, pos, name)
M = size(pos,1);
if isnumeric(rho) && isreal(rho) && isscalar(rho) && isfinite(rho)
    q = repmat(rho,M,1);
elseif isa(rho,'function_handle')
    q = rho(pos);
    if isscalar(q), q = repmat(q,M,1); end
    if ~(isnumeric(q) && isreal(q) && numel(q) == M && all(isfinite(q(:))))
        error('em:src:volCharge:InvalidDensityOutput', ...
            '%s(pos) must return one finite real value per position.', name);
    end
    q = q(:);
else
    error('em:src:volCharge:InvalidDensity', ...
        '%s must be a finite real scalar or a function handle.', name);
end
end
