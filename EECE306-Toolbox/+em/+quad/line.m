function [pos, w, tangent] = line(curve, tspan, N, rule)
%LINE Discretize a parameterized curve with arc-length weights.
%   [POS, W] = EM.QUAD.LINE(CURVE, TSPAN, N, RULE) returns Nx3 Cartesian
%   positions and Nx1 weights in meters. CURVE maps a scalar parameter to
%   a 1x3 position. The default rule is 'midpoint'.
%   Optional third output TANGENT is Nx3 dr/dt before normalization or
%   quadrature weighting. Existing two-output callers are unchanged.
%
%   Example:
%       c = @(t) [cos(t) sin(t) 0];
%       [pos,w] = em.quad.line(c, [0 2*pi], 200);
%
%   See also EM.QUAD.NODES, EM.QUAD.SURF, EM.QUAD.VOL.

if nargin < 4 || isempty(rule), rule = 'midpoint'; end
validateParametrization(curve, 'curve');
validateSpan(tspan, 'tspan');
[t, wnode] = em.quad.nodes(tspan(1), tspan(2), N, rule);
dt = 1e-6 * (tspan(2)-tspan(1));
pos = evaluateCurve(curve, t);
dpos = (evaluateCurve(curve, t+dt) - evaluateCurve(curve, t-dt))/(2*dt);
tangent = dpos;
w = sqrt(sum(dpos.^2,2)) .* wnode;
end

function pos = evaluateCurve(curve, t)
M = numel(t);
try
    candidate = curve(t);
catch
    candidate = [];
end
if isnumeric(candidate) && isreal(candidate) && isequal(size(candidate), [M 3])
    pos = candidate;
else
    pos = zeros(M,3);
    for k = 1:M
        value = curve(t(k));
        if ~(isnumeric(value) && isreal(value) && numel(value) == 3)
            error('em:quad:line:InvalidCurveOutput', ...
                'curve must return one finite Cartesian 1x3 position.');
        end
        pos(k,:) = reshape(value,1,3);
    end
end
if any(~isfinite(pos(:)))
    error('em:quad:line:InvalidCurveOutput', ...
        'curve must return finite Cartesian positions.');
end
end

function validateParametrization(f, name)
if ~isa(f, 'function_handle')
    error('em:quad:line:InvalidCurve', '%s must be a function handle.', name);
end
end

function validateSpan(span, name)
if ~(isnumeric(span) && isreal(span) && isequal(size(span), [1 2]) ...
        && all(isfinite(span)) && span(2) > span(1))
    error('em:quad:line:InvalidSpan', ...
        '%s must be a finite increasing 1x2 array.', name);
end
end
