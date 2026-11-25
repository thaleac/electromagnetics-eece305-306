function g = grad(f, r, h)
%GRAD Central-difference gradient of a scalar field.
%   G = EM.OP.GRAD(F,R,H) returns Nx3 at Nx3 positions R (m). F must
%   accept all positions together and return finite real Nx1 values.
%   H defaults to 1e-5 m. Each component uses (f(r+h*ek)-f(r-h*ek))/(2*h).
%   Units are scalar-field units per meter. Smooth-field truncation is
%   O(h^2); excessively small steps amplify floating-point cancellation.
%
%   Example:
%       g = em.op.grad(@(r) sum(r.^2,2),[1 2 3]); % [2 4 6]
%
%   See also EM.FIELD.V, EM.OP.DIV.
if nargin<3 || isempty(h), h=1e-5; end
if ~isa(f,'function_handle'), error('em:op:grad:InvalidField','f must be a function handle.'); end
if ~(isnumeric(r) && isreal(r) && ismatrix(r) && size(r,2)==3 && all(isfinite(r(:))))
    error('em:op:grad:InvalidPosition','r must be a finite real Nx3 position array.');
end
if ~(isnumeric(h) && isreal(h) && isscalar(h) && isfinite(h) && h>0)
    error('em:op:grad:InvalidStep','h must be a finite positive scalar in meters.');
end
g=zeros(size(r));
if isempty(r), return; end
for k=1:3
    rp=r; rm=r; rp(:,k)=rp(:,k)+h; rm(:,k)=rm(:,k)-h;
    if any(rp(:,k)==rm(:,k)), error('em:op:grad:UnresolvedStep','h cannot perturb r in direction %d.',k); end
    fp=f(rp); fm=f(rm);
    for value={fp,fm}
        a=value{1};
        if ~(isnumeric(a) && isreal(a) && isequal(size(a),[size(r,1) 1]) && all(isfinite(a(:))))
            error('em:op:grad:InvalidFieldOutput','f must return finite real Nx1 values.');
        end
    end
    g(:,k)=(fp-fm)/(2*h);
end
end
