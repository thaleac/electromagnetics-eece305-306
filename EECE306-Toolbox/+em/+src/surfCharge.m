function s = surfCharge(rhoS, surface, uspan, vspan, Nu, Nv, rule)
%SURFCHARGE Construct a discretized surface-charge source.
%   S = EM.SRC.SURFCHARGE(RHOS, SURFACE, USPAN, VSPAN, NU, NV, RULE)
%   accepts a scalar density in C/m^2 or a function handle RHOS(POS). The
%   default quadrature rule is 'midpoint'.
%
%   Example:
%       d = @(u,v) [u.*cos(v) u.*sin(v) zeros(size(u))];
%       s = em.src.surfCharge(1e-9,d,[0 1],[0 2*pi],40,80);
%
%   See also EM.QUAD.SURF, EM.SRC.LINECHARGE, EM.SRC.VOLCHARGE.

if nargin < 7 || isempty(rule), rule = 'midpoint'; end
[pos,w] = em.quad.surf(surface, uspan, vspan, Nu, Nv, rule);
q = evaluateDensity(rhoS, pos, 'rhoS');
s = struct('type','charge', 'pos',pos, 'w',w, 'q',q, ...
    'label','surface charge');
end

function q = evaluateDensity(rho, pos, name)
M = size(pos,1);
if isnumeric(rho) && isreal(rho) && isscalar(rho) && isfinite(rho)
    q = repmat(rho,M,1);
elseif isa(rho,'function_handle')
    q = rho(pos);
    if isscalar(q), q = repmat(q,M,1); end
    if ~(isnumeric(q) && isreal(q) && numel(q) == M && all(isfinite(q(:))))
        error('em:src:surfCharge:InvalidDensityOutput', ...
            '%s(pos) must return one finite real value per position.', name);
    end
    q = q(:);
else
    error('em:src:surfCharge:InvalidDensity', ...
        '%s must be a finite real scalar or a function handle.', name);
end
end
