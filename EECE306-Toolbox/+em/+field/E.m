function Ev = E(s, r)
%E Evaluate the electric field of a discretized charge source.
%   EV = EM.FIELD.E(S, R) evaluates source S at the Nx3 Cartesian
%   positions R and returns the Nx3 electric field in V/m.
%
%   Example:
%       s = em.src.pointCharge(1e-9, [0 0 0]);
%       Ev = em.field.E(s, [1 0 0]);
%
%   See also EM.SRC.POINTCHARGE, EM.SRC.MERGE, EM.CONST.EPS0.

validateChargeSource(s);
if ~(isnumeric(r) && isreal(r) && ismatrix(r) && size(r, 2) == 3 ...
        && all(isfinite(r(:))))
    error('em:field:E:InvalidPosition', ...
        'r must be a finite real Nx3 Cartesian position array.');
end

N = size(r, 1);
M = size(s.pos, 1);
Ev = zeros(N, 3);
if N == 0 || M == 0
    return
end

% Limit temporary displacement arrays to roughly one million pairs.
maxPairsPerChunk = 1e6;
chunkSize = max(1, floor(maxPairsPerChunk / M));
qw = (s.q(:) .* s.w(:)).';
scale = max([1; abs(r(:)); abs(s.pos(:))]);
singularTol = 1e-12 * scale;

for first = 1:chunkSize:N
    last = min(N, first + chunkSize - 1);
    rc = r(first:last, :);
    nc = size(rc, 1);
    d = reshape(rc, [nc 1 3]) - reshape(s.pos, [1 M 3]);
    R2 = sum(d.^2, 3);
    [obsLocal, elemIndex] = find(R2 <= singularTol^2, 1);
    if ~isempty(elemIndex)
        error('em:field:E:Singularity', ...
            ['Observation point %d lies at or too near source element ' ...
             '%d; the electric field is singular.'], ...
            first + obsLocal - 1, elemIndex);
    end
    kernel = qw ./ (R2 .* sqrt(R2));
    block = sum(d .* reshape(kernel, [nc M 1]), 2);
    Ev(first:last, :) = reshape(block, [nc 3]);
end

Ev = Ev / (4*pi*em.const.eps0());
end

function validateChargeSource(s)
if ~isstruct(s) || ~isfield(s, 'type') || ~strcmp(s.type, 'charge')
    error('em:field:E:InvalidSource', ...
        's must be a charge source struct.');
end
required = {'pos', 'w', 'q'};
for k = 1:numel(required)
    if ~isfield(s, required{k})
        error('em:field:E:InvalidSource', ...
            's is missing required field %s.', required{k});
    end
end
M = size(s.pos, 1);
if ~isnumeric(s.pos) || size(s.pos, 2) ~= 3 ...
        || ~isnumeric(s.w) || numel(s.w) ~= M ...
        || ~isnumeric(s.q) || numel(s.q) ~= M
    error('em:field:E:InvalidSource', ...
        's.pos, s.w, and s.q must have compatible M-element shapes.');
end
end
