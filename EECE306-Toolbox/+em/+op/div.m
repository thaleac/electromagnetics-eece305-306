function d = div(F, r, h)
%DIV Central-difference divergence of a vector field.
%   D = EM.OP.DIV(F, R, H) returns Nx1 divergence values at Nx3
%   Cartesian positions R (m). F accepts Nx3 and returns Nx3. The
%   positive scalar step H defaults to 1e-5 m. For each direction k,
%   only component k is differentiated: [Fk(r+h*ek)-Fk(r-h*ek)]/(2*h).
%   Output units are field units per meter (C/m^3 for D).
%
%   Truncation error is O(h^2) for smooth fields. A discretized volume
%   source consists of point elements: h much smaller than their spacing
%   resolves charge-free gaps rather than the continuum density. Use a
%   refinement and step-size study before interpreting interior values.
%
%   Example:
%       F = @(r) r;
%       d = em.op.div(F, [1 2 3]); % approximately 3
%
%   See also EM.FIELD.D, EM.FIELD.FLUX.

if nargin < 3 || isempty(h), h = 1e-5; end
if ~isa(F,'function_handle')
    error('em:op:div:InvalidField','F must be a function handle.');
end
if ~(isnumeric(r) && isreal(r) && ismatrix(r) && size(r,2) == 3 ...
        && all(isfinite(r(:))))
    error('em:op:div:InvalidPosition','r must be a finite real Nx3 position array.');
end
if ~(isnumeric(h) && isreal(h) && isscalar(h) && isfinite(h) && h > 0)
    error('em:op:div:InvalidStep','h must be a finite positive real scalar in meters.');
end
d = zeros(size(r,1),1);
if isempty(r), return; end
for k = 1:3
    rp = r; rm = r;
    rp(:,k) = rp(:,k) + h;
    rm(:,k) = rm(:,k) - h;
    if any(rp(:,k) == rm(:,k))
        error('em:op:div:UnresolvedStep','h is too small to perturb r in direction %d.',k);
    end
    fp = F(rp); fm = F(rm);
    validateField(fp,size(r),'F(r+h*ek)');
    validateField(fm,size(r),'F(r-h*ek)');
    d = d + (fp(:,k)-fm(:,k))/(2*h);
end
end

function validateField(value,shape,name)
if ~(isnumeric(value) && isreal(value) && isequal(size(value),shape) ...
        && all(isfinite(value(:))))
    error('em:op:div:InvalidFieldOutput', ...
        '%s must return a finite real Nx3 vector field.',name);
end
end
