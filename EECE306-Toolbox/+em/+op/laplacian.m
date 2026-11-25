function L = laplacian(f, r, h)
%LAPLACIAN Central-difference Laplacian of a scalar field.
%   L = EM.OP.LAPLACIAN(F,R,H) returns Nx1 values at the Nx3 Cartesian
%   positions R. F must accept all positions together and return finite
%   real Nx1 values. H defaults to 1e-5 m. The operator sums the second
%   central difference in each Cartesian direction.
%
%   Example:
%       L = em.op.laplacian(@(r) sum(r.^2,2), [1 2 3]); % approximately 6
%
%   See also EM.OP.GRAD, EM.OP.DIV.

if nargin < 3 || isempty(h), h = 1e-5; end
if ~isa(f, 'function_handle')
    error('em:op:laplacian:InvalidField', 'f must be a function handle.');
end
if ~(isnumeric(r) && isreal(r) && ismatrix(r) && size(r,2) == 3 ...
        && all(isfinite(r(:))))
    error('em:op:laplacian:InvalidPosition', ...
        'r must be a finite real Nx3 position array.');
end
if ~(isnumeric(h) && isreal(h) && isscalar(h) && isfinite(h) && h > 0)
    error('em:op:laplacian:InvalidStep', ...
        'h must be a finite positive scalar in meters.');
end

L = zeros(size(r,1), 1);
if isempty(r), return; end
f0 = f(r);
validateField(f0, size(r,1), 'f(r)');
for k = 1:3
    rp = r; rm = r;
    rp(:,k) = rp(:,k) + h;
    rm(:,k) = rm(:,k) - h;
    if any(rp(:,k) == rm(:,k))
        error('em:op:laplacian:UnresolvedStep', ...
            'h cannot perturb r in direction %d.', k);
    end
    fp = f(rp); fm = f(rm);
    validateField(fp, size(r,1), 'f(r+h*ek)');
    validateField(fm, size(r,1), 'f(r-h*ek)');
    L = L + (fp - 2*f0 + fm) / h^2;
end
end

function validateField(value, N, name)
if ~(isnumeric(value) && isreal(value) && isequal(size(value), [N 1]) ...
        && all(isfinite(value(:))))
    error('em:op:laplacian:InvalidFieldOutput', ...
        '%s must return finite real Nx1 values.', name);
end
end
