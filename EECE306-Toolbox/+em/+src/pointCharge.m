function s = pointCharge(Q, r0)
%POINTCHARGE Construct a point-charge source.
%   S = EM.SRC.POINTCHARGE(Q, R0) returns a charge-source struct for a
%   charge Q in coulombs at the 1x3 Cartesian position R0 in meters.
%
%   Example:
%       s = em.src.pointCharge(1e-9, [0 0 0]);
%
%   See also EM.SRC.MERGE, EM.FIELD.E.

if ~(isnumeric(Q) && isreal(Q) && isscalar(Q) && isfinite(Q))
    error('em:src:pointCharge:InvalidQ', ...
        'Q must be a finite real scalar in coulombs.');
end
if ~(isnumeric(r0) && isreal(r0) && isequal(size(r0), [1 3]) ...
        && all(isfinite(r0)))
    error('em:src:pointCharge:InvalidPosition', ...
        'r0 must be a finite real 1x3 Cartesian position.');
end

s = struct('type', 'charge', 'pos', r0, 'w', 1, 'q', Q, ...
    'label', 'point charge');
end
