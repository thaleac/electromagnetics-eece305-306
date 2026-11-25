function th = angle(A, B)
%ANGLE Interior angle between two sets of vectors
%   TH = EM.VEC.ANGLE(A, B) returns the Nx1 interior angle in radians
%   between the Nx3 arrays A and B
%
%   Example:
%       em.vec.angle([1 0 0], [0 1 0]) % returns 1.5708 (pi/2)
%
%   See also EM.VEC.MAG, EM.VEC.FROMTO

    if ~ismatrix(A) || size(A,2) ~= 3
        error('First input must be an Nx3 array');
    end

    if ~ismatrix(B) || size(B,2) ~= 3
        error('Second input must be an Nx3 array');
    end

    if size(A,1) ~= size(B,1)
        error('Second input must have the same number of rows as first input');
    end

    dp = sum(A .* B, 2);
    m = em.vec.mag(A) .* em.vec.mag(B);
    x = max(-1, min(1, dp ./ m));
    th = acos(x);
end
