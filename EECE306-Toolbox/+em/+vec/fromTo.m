function d = fromTo(P, Q)
%FROMTO Displacement vector between two points
%   D = EM.VEC.FROMTO(P, Q) returns the Nx3 displacement vector 
%   from point P to point Q
%
%   Example:
%       em.vec.fromTo([1 1 1], [2 3 4]) % returns [1 2 3]
%
%   See also EM.VEC.MAG, EM.VEC.UNIT

    if ~ismatrix(P) || size(P,2) ~= 3
        error('First input must be an Nx3 array');
    end

    if ~ismatrix(Q) || size(Q,2) ~= 3
        error('Second input must be an Nx3 array');
    end

    if size(P,1) ~= size(Q,1)
        error('Second input must have the same number of rows as first input');
    end

    d = Q - P;
end
