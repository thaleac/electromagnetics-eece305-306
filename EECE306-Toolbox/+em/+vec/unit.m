function u = unit(A)
%UNIT Unit vector of each row of A
%   U = EM.VEC.UNIT(A) returns the Nx3 array of unit vectors of the 
%   Nx3 input A. A zero row produces a row of NaN
%
%   Example:
%       em.vec.unit([3 4 0]) % returns [0.6 0.8 0]
%
%   See also EM.VEC.MAG

    if ~ismatrix(A) || size(A,2) ~= 3
        error('Input must be an Nx3 array');
    end

    u = A ./ em.vec.mag(A);
end
