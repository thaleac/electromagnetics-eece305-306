function m = mag(A)
% MAG(A) - The magnitude of a vector
%       M = MAG(A) returns the Nx1 magnitude array for
%       the Nx3 input array of vector A
%
%   Example:
%       m = em.vec.mag([3,4,0]) % returns 5
%   See also EM.VEC.UNIT.

    if ~ismatrix(A) || size(A,2) ~= 3
        error('A must be an Nx3 array.');
    end

    m = sqrt(sum(A.^2, 2));
end
