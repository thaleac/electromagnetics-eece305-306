function s = lineCharge(rhoL, curve, tspan, N, rule)
%LINECHARGE Construct a discretized line-charge source.
%   S = EM.SRC.LINECHARGE(RHOL, CURVE, TSPAN, N, RULE) constructs a charge
%   source from a scalar line density RHOL in C/m or a function handle
%   RHOL(POS). The default quadrature rule is 'midpoint'.
%
%   Example:
%       c = @(t) [0 0 t];
%       s = em.src.lineCharge(1e-9, c, [-1 1], 200);
%
%   See also EM.QUAD.LINE, EM.SRC.SURFCHARGE, EM.SRC.VOLCHARGE.

if nargin < 5 || isempty(rule), rule = 'midpoint'; end
[pos,w] = em.quad.line(curve, tspan, N, rule);
q = evaluateDensity(rhoL, pos, 'rhoL');
s = struct('type','charge', 'pos',pos, 'w',w, 'q',q, ...
    'label','line charge');
end

function q = evaluateDensity(rho, pos, name)
M = size(pos,1);
if isnumeric(rho) && isreal(rho) && isscalar(rho) && isfinite(rho)
    q = repmat(rho,M,1);
elseif isa(rho,'function_handle')
    q = rho(pos);
    if isscalar(q), q = repmat(q,M,1); end
    if ~(isnumeric(q) && isreal(q) && numel(q) == M && all(isfinite(q(:))))
        error('em:src:lineCharge:InvalidDensityOutput', ...
            '%s(pos) must return one finite real value per position.', name);
    end
    q = q(:);
else
    error('em:src:lineCharge:InvalidDensity', ...
        '%s must be a finite real scalar or a function handle.', name);
end
end
