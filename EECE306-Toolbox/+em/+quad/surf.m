function [pos, w] = surf(surface, uspan, vspan, Nu, Nv, rule)
%SURF Discretize a parameterized surface with area weights.
%   [POS, W] = EM.QUAD.SURF(SURFACE, USPAN, VSPAN, NU, NV, RULE) returns
%   Mx3 Cartesian positions and Mx1 weights in square meters, M=NU*NV.
%   SURFACE maps scalar parameters (u,v) to a 1x3 position. The default
%   rule is 'midpoint'.
%
%   Example:
%       s = @(u,v) [u v 0];
%       [pos,w] = em.quad.surf(s, [0 1], [0 1], 20, 20);
%
%   See also EM.QUAD.NODES, EM.QUAD.LINE, EM.QUAD.VOL.

if nargin < 6 || isempty(rule), rule = 'midpoint'; end
if ~isa(surface, 'function_handle')
    error('em:quad:surf:InvalidSurface', 'surface must be a function handle.');
end
validateSpan(uspan, 'uspan');
validateSpan(vspan, 'vspan');
[u, wu] = em.quad.nodes(uspan(1), uspan(2), Nu, rule);
[v, wv] = em.quad.nodes(vspan(1), vspan(2), Nv, rule);
[U,V] = ndgrid(u,v);
uf = U(:); vf = V(:);
du = 1e-6 * (uspan(2)-uspan(1));
dv = 1e-6 * (vspan(2)-vspan(1));
pos = evaluateSurface(surface, uf, vf);
ru = (evaluateSurface(surface, uf+du, vf) ...
    - evaluateSurface(surface, uf-du, vf))/(2*du);
rv = (evaluateSurface(surface, uf, vf+dv) ...
    - evaluateSurface(surface, uf, vf-dv))/(2*dv);
jac = sqrt(sum(cross(ru,rv,2).^2,2));
[WU,WV] = ndgrid(wu,wv);
w = jac .* (WU(:).*WV(:));
end

function pos = evaluateSurface(surface, u, v)
M = numel(u);
try
    candidate = surface(u,v);
catch
    candidate = [];
end
if isnumeric(candidate) && isreal(candidate) && isequal(size(candidate), [M 3])
    pos = candidate;
else
    pos = zeros(M,3);
    for k = 1:M
        value = surface(u(k),v(k));
        if ~(isnumeric(value) && isreal(value) && numel(value) == 3)
            error('em:quad:surf:InvalidSurfaceOutput', ...
                'surface must return one finite Cartesian 1x3 position.');
        end
        pos(k,:) = reshape(value,1,3);
    end
end
if any(~isfinite(pos(:)))
    error('em:quad:surf:InvalidSurfaceOutput', ...
        'surface must return finite Cartesian positions.');
end
end

function validateSpan(span, name)
if ~(isnumeric(span) && isreal(span) && isequal(size(span), [1 2]) ...
        && all(isfinite(span)) && span(2) > span(1))
    error('em:quad:surf:InvalidSpan', ...
        '%s must be a finite increasing 1x2 array.', name);
end
end
