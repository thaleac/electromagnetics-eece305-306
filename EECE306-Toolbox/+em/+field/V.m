function potential = V(s, r)
%V Electric potential of a discretized finite charge source.
%   POTENTIAL = EM.FIELD.V(S,R) returns Nx1 volts at Nx3 positions R (m),
%   using sum(q.*w/R)/(4*pi*eps0) and the reference V=0 at infinity.
%   Observation points within 1e-12 times the coordinate scale of an
%   element raise a singularity error naming that element.
%
%   The required element list cannot reveal a family's infinite-extent
%   intent. Optional S.extentGrowsWithN (logical scalar, default false)
%   declares increasing physical extent toward an infinite distribution.
%   If true, a warning accompanies the finite truncation's returned V;
%   that value does not establish an infinite source's reference at infinity.
%   Fixed-geometry refinement, large element counts, and large potentials
%   alone do not trigger this warning. Set this flag after construction or
%   merging: existing constructors do not propagate optional metadata.
%
%   Example:
%       s = em.src.pointCharge(1e-9,[0 0 0]);
%       potential = em.field.V(s,[1 0 0]);
%
%   See also EM.FIELD.E, EM.OP.GRAD, EM.SRC.LINECHARGE.
if ~(isstruct(s) && isscalar(s) && isfield(s,'type') && strcmp(s.type,'charge'))
    error('em:field:V:InvalidSource','s must be a scalar charge source struct.');
end
names = {'pos','q','w'};
for k = 1:3
    if ~isfield(s,names{k})
        error('em:field:V:InvalidSource','s is missing %s.',names{k});
    end
    a = s.(names{k});
    if ~(isnumeric(a) && isreal(a) && all(isfinite(a(:))))
        error('em:field:V:InvalidSource','s.%s must be finite real numeric data.',names{k});
    end
end
M = size(s.pos,1);
if ~ismatrix(s.pos) || size(s.pos,2) ~= 3 || numel(s.q) ~= M || numel(s.w) ~= M
    error('em:field:V:InvalidSource','s.pos must be Mx3; s.q and s.w must have M entries.');
end
if ~(isnumeric(r) && isreal(r) && ismatrix(r) && size(r,2)==3 && all(isfinite(r(:))))
    error('em:field:V:InvalidPosition','r must be a finite real Nx3 position array.');
end
if isfield(s,'extentGrowsWithN')
    if ~(islogical(s.extentGrowsWithN) && isscalar(s.extentGrowsWithN))
        error('em:field:V:InvalidExtent','s.extentGrowsWithN must be a logical scalar.');
    end
    if s.extentGrowsWithN
        warning('em:field:V:InfiniteReference', ...
            's declares extent growing with element count; V=0 at infinity may not exist. Returned V describes only this finite truncation.');
    end
end
potential = zeros(size(r,1),1);
if M==0 || isempty(r), return; end
qw = (s.q(:).*s.w(:)).';
tol = 1e-12*max([1;abs(r(:));abs(s.pos(:))]);
chunk = max(1,floor(1e6/M));
for first = 1:chunk:size(r,1)
    last = min(size(r,1),first+chunk-1);
    d = reshape(r(first:last,:),[],1,3)-reshape(s.pos,1,M,3);
    R2 = sum(d.^2,3);
    [obs,elem] = find(R2 <= tol^2,1);
    if ~isempty(elem)
        error('em:field:V:Singularity','Observation point %d lies at or too near source element %d; V is singular.',first+obs-1,elem);
    end
    potential(first:last) = sum(qw./sqrt(R2),2)/(4*pi*em.const.eps0());
end
end
