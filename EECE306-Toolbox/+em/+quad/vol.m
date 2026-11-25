function [pos, w] = vol(volume, uspan, vspan, wspan, Nu, Nv, Nw, rule)
%VOL Discretize a parameterized volume with volume weights.
%   [POS, W] = EM.QUAD.VOL(VOLUME, USPAN, VSPAN, WSPAN, NU, NV, NW,
%   RULE) returns Mx3 Cartesian positions and Mx1 weights in cubic meters,
%   M=NU*NV*NW. The default rule is 'midpoint'.
%
%   Example:
%       v = @(u,v,w) [u v w];
%       [pos,wt] = em.quad.vol(v,[0 1],[0 1],[0 1],10,10,10);
%
%   See also EM.QUAD.NODES, EM.QUAD.LINE, EM.QUAD.SURF.

if nargin < 8 || isempty(rule), rule = 'midpoint'; end
if ~isa(volume, 'function_handle')
    error('em:quad:vol:InvalidVolume', 'volume must be a function handle.');
end
validateSpan(uspan, 'uspan');
validateSpan(vspan, 'vspan');
validateSpan(wspan, 'wspan');
[u, wu] = em.quad.nodes(uspan(1), uspan(2), Nu, rule);
[v, wv] = em.quad.nodes(vspan(1), vspan(2), Nv, rule);
[q, ww] = em.quad.nodes(wspan(1), wspan(2), Nw, rule);
[U,V,Q] = ndgrid(u,v,q);
uf = U(:); vf = V(:); qf = Q(:);
du = 1e-6 * (uspan(2)-uspan(1));
dv = 1e-6 * (vspan(2)-vspan(1));
dq = 1e-6 * (wspan(2)-wspan(1));
pos = evaluateVolume(volume, uf, vf, qf);
ru = (evaluateVolume(volume, uf+du, vf, qf) ...
    - evaluateVolume(volume, uf-du, vf, qf))/(2*du);
rv = (evaluateVolume(volume, uf, vf+dv, qf) ...
    - evaluateVolume(volume, uf, vf-dv, qf))/(2*dv);
rw = (evaluateVolume(volume, uf, vf, qf+dq) ...
    - evaluateVolume(volume, uf, vf, qf-dq))/(2*dq);
jac = abs(dot(ru, cross(rv,rw,2), 2));
[WU,WV,WW] = ndgrid(wu,wv,ww);
w = jac .* (WU(:).*WV(:).*WW(:));
end

function pos = evaluateVolume(volume, u, v, q)
M = numel(u);
try
    candidate = volume(u,v,q);
catch
    candidate = [];
end
if isnumeric(candidate) && isreal(candidate) && isequal(size(candidate), [M 3])
    pos = candidate;
else
    pos = zeros(M,3);
    for k = 1:M
        value = volume(u(k),v(k),q(k));
        if ~(isnumeric(value) && isreal(value) && numel(value) == 3)
            error('em:quad:vol:InvalidVolumeOutput', ...
                'volume must return one finite Cartesian 1x3 position.');
        end
        pos(k,:) = reshape(value,1,3);
    end
end
if any(~isfinite(pos(:)))
    error('em:quad:vol:InvalidVolumeOutput', ...
        'volume must return finite Cartesian positions.');
end
end

function validateSpan(span, name)
if ~(isnumeric(span) && isreal(span) && isequal(size(span), [1 2]) ...
        && all(isfinite(span)) && span(2) > span(1))
    error('em:quad:vol:InvalidSpan', ...
        '%s must be a finite increasing 1x2 array.', name);
end
end
