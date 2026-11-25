function c = curl(F, r, h)
%CURL Central-difference curl of a Cartesian vector field.
%   C = EM.OP.CURL(F,R,H) returns Nx3 values at Nx3 positions R (m).
%   F accepts Nx3 positions and returns finite real Nx3 vectors.
%   H is a positive scalar step in meters and defaults to 1e-5.
%   Smooth fields have O(H^2) truncation error.
%
%   Example:
%       c = em.op.curl(@(r) [-r(:,2), r(:,1), 0*r(:,3)], [1 2 3]);
%       % c is approximately [0 0 2], in field units per meter.
%
%   See also EM.OP.DIV, EM.OP.GRAD.
if nargin < 3 || isempty(h), h = 1e-5; end
if ~isa(F,'function_handle')
    error('em:op:curl:InvalidField','F must be a function handle.');
end
if ~(isnumeric(r) && isreal(r) && ismatrix(r) && size(r,2) == 3 && all(isfinite(r(:))))
    error('em:op:curl:InvalidPosition','r must be a finite real Nx3 array.');
end
if ~(isnumeric(h) && isreal(h) && isscalar(h) && isfinite(h) && h > 0)
    error('em:op:curl:InvalidStep','h must be a finite positive real scalar.');
end
c = zeros(size(r));
if isempty(r), return; end

derivative = cell(1,3);
for k = 1:3
    rp = double(r); rm = rp;
    rp(:,k) = rp(:,k) + h;
    rm(:,k) = rm(:,k) - h;
    if any(rp(:,k) == rm(:,k))
        error('em:op:curl:UnresolvedStep','h cannot perturb r in direction %d.',k);
    end
    fp = F(rp); fm = F(rm);
    for output = {fp,fm}
        value = output{1};
        if ~(isnumeric(value) && isreal(value) && isequal(size(value),size(r)) ...
                && all(isfinite(value(:))))
            error('em:op:curl:InvalidFieldOutput','F must return finite real Nx3 vectors.');
        end
    end
    derivative{k} = (double(fp)-double(fm))/(2*h);
end
dx = derivative{1}; dy = derivative{2}; dz = derivative{3};
c = [dy(:,3)-dz(:,2), dz(:,1)-dx(:,3), dx(:,2)-dy(:,1)];
end
