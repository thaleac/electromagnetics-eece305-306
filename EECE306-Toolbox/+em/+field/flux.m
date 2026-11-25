function Phi = flux(F, surface, uspan, vspan, Nu, Nv, rule)
%FLUX Integrate a vector field over an oriented parameterized surface.
%   PHI = EM.FIELD.FLUX(F, SURFACE, USPAN, VSPAN, NU, NV, RULE) returns
%   the scalar integral of F dot dS. F accepts Mx3 Cartesian positions
%   and returns Mx3 vectors. SURFACE maps scalar (u,v) to a 1x3 position
%   in meters; vectorized Mx3 output is also supported. RULE defaults to
%   'midpoint'; the four EM.QUAD.NODES rules are supported.
%
%   Orientation is ru cross rv, in parameter order (u,v). Swapping the
%   parameters reverses the flux. Outward orientation and closure are
%   the caller's responsibility. Scalar area weights cannot encode this
%   direction. Tangents use central differences with steps 1e-6 times
%   the respective parameter spans. The surface must extend smoothly to
%   these nearby points. Units are field units times m^2 (C for D).
%   Source singularities, including the element index, propagate from F.
%
%   Example:
%       s = em.src.pointCharge(1e-9, [0 0 0]);
%       sph = @(th,ph) [sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
%       Phi = em.field.flux(@(r) em.field.D(s,r), sph, ...
%           [0 pi], [0 2*pi], 32, 64, 'gauss');
%
%   See also EM.QUAD.NODES, EM.QUAD.SURF, EM.FIELD.D.

if nargin < 7 || isempty(rule), rule = 'midpoint'; end
if ~isa(F, 'function_handle')
    error('em:field:flux:InvalidField', 'F must be a function handle.');
end
if ~isa(surface, 'function_handle')
    error('em:field:flux:InvalidSurface', 'surface must be a function handle.');
end
validateSpan(uspan, 'uspan');
validateSpan(vspan, 'vspan');
validateCount(Nu, 'Nu');
validateCount(Nv, 'Nv');
[u,wu] = em.quad.nodes(uspan(1),uspan(2),Nu,rule);
[v,wv] = em.quad.nodes(vspan(1),vspan(2),Nv,rule);
[U,V] = ndgrid(u,v);
u = U(:); v = V(:);
du = 1e-6 * diff(uspan); dv = 1e-6 * diff(vspan);
pos = evaluateSurface(surface,u,v);
ru = (evaluateSurface(surface,u+du,v) ...
    - evaluateSurface(surface,u-du,v))/(2*du);
rv = (evaluateSurface(surface,u,v+dv) ...
    - evaluateSurface(surface,u,v-dv))/(2*dv);
areaVector = cross(ru,rv,2);
values = F(pos);
if ~(isnumeric(values) && isreal(values) ...
        && isequal(size(values),size(pos)) && all(isfinite(values(:))))
    error('em:field:flux:InvalidFieldOutput', ...
        'F(pos) must return a finite real Mx3 vector field.');
end
[WU,WV] = ndgrid(wu,wv);
Phi = sum(sum(values.*areaVector,2) .* WU(:) .* WV(:));
end

function pos = evaluateSurface(surface,u,v)
M = numel(u);
try
    candidate = surface(u,v);
catch
    candidate = [];
end
if isnumeric(candidate) && isreal(candidate) && isequal(size(candidate),[M 3])
    pos = candidate;
else
    pos = zeros(M,3);
    for k = 1:M
        value = surface(u(k),v(k));
        if ~(isnumeric(value) && isreal(value) && isequal(size(value),[1 3]))
            error('em:field:flux:InvalidSurfaceOutput', ...
                'surface must return a finite real 1x3 position.');
        end
        pos(k,:) = value;
    end
end
if any(~isfinite(pos(:)))
    error('em:field:flux:InvalidSurfaceOutput', ...
        'surface must return finite Cartesian positions.');
end
end

function validateSpan(span,name)
if ~(isnumeric(span) && isreal(span) && isequal(size(span),[1 2]) ...
        && all(isfinite(span)) && span(2) > span(1))
    error('em:field:flux:InvalidSpan', ...
        '%s must be a finite increasing 1x2 array.',name);
end
end

function validateCount(N,name)
if ~(isnumeric(N) && isreal(N) && isscalar(N) && isfinite(N) ...
        && N == floor(N) && N >= 1)
    error('em:field:flux:InvalidCount','%s must be a positive integer scalar.',name);
end
end
