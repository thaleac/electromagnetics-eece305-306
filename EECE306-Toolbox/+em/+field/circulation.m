function C = circulation(F, curve, tspan, N, rule)
%CIRCULATION Integrate a vector field along a directed open or closed path.
%   C = EM.FIELD.CIRCULATION(F,CURVE,TSPAN,N,RULE) returns scalar
%   sum(F(r(t)).*rprime(t)*quadratureWeight). F accepts Nx3 positions and
%   returns Nx3 vectors. CURVE accepts scalar or Nx1 t and returns 1x3 or
%   Nx3 positions in meters. RULE defaults to 'midpoint'; all four rules
%   of EM.QUAD.NODES are supported. Units are field units times meters.
%
%   Tangents retain direction and use a central parameter step of 1e-6
%   times diff(TSPAN). The curve must extend to nearby parameters. At a
%   corner this gives the mean of the one-sided tangents. For piecewise
%   smooth paths, split at corners or align midpoint panels with them;
%   quadrature spanning unresolved corners requires refinement. Closure
%   is not imposed. Reverse CURVE's parameterization to reverse direction.
%   Source-singularity errors (including element indices) propagate from F.
%
%   Example:
%       c = @(t) [2*cos(t) 2*sin(t) zeros(size(t))];
%       C = em.field.circulation(@(r) [-r(:,2) r(:,1) 0*r(:,1)], ...
%           c,[0 2*pi],400); % approximately 8*pi
%
%   See also EM.QUAD.NODES, EM.FIELD.FLUX, EM.FIELD.V.
if nargin<5 || isempty(rule), rule='midpoint'; end
if ~isa(F,'function_handle'), error('em:field:circulation:InvalidField','F must be a function handle.'); end
if ~isa(curve,'function_handle'), error('em:field:circulation:InvalidCurve','curve must be a function handle.'); end
if ~(isnumeric(tspan) && isreal(tspan) && isequal(size(tspan),[1 2]) && all(isfinite(tspan)) && tspan(2)>tspan(1))
    error('em:field:circulation:InvalidSpan','tspan must be a finite increasing 1x2 array.');
end
[t,w]=em.quad.nodes(tspan(1),tspan(2),N,rule);
dt=1e-6*diff(tspan);
pos=evaluateCurve(curve,t);
tangent=(evaluateCurve(curve,t+dt)-evaluateCurve(curve,t-dt))/(2*dt);
values=F(pos);
if ~(isnumeric(values) && isreal(values) && isequal(size(values),size(pos)) && all(isfinite(values(:))))
    error('em:field:circulation:InvalidFieldOutput','F must return finite real Nx3 vectors.');
end
C=sum(sum(values.*tangent,2).*w);
end

function pos=evaluateCurve(curve,t)
try
    pos=curve(t);
catch
    pos=[];
end
if ~(isnumeric(pos) && isreal(pos) && isequal(size(pos),[numel(t) 3]))
    pos=zeros(numel(t),3);
    for k=1:numel(t)
        a=curve(t(k));
        if ~(isnumeric(a) && isreal(a) && isequal(size(a),[1 3]))
            error('em:field:circulation:InvalidCurveOutput','curve must return finite real 1x3 or Nx3 positions.');
        end
        pos(k,:)=a;
    end
end
if any(~isfinite(pos(:)))
    error('em:field:circulation:InvalidCurveOutput','curve must return finite real positions.');
end
end
