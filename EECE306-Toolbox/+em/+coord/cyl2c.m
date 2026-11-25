function r = cyl2c(p)
%CYL2C - From cylindrical to cartesian coordinates
%   R = EM.COORD.CYL2C(P) converts the Nx3 cylindrical points
%   P = [rho phi z] into Cartesian coordinates R = [x y z]. Angles are
%   in radians
%
%   Example:
%       r = em.coord.cyl2c([1 pi/4 0]) % returns [sqrt(2)/2  sqrt(2)/2  0]
%
%   See also EM.COORD.C2CYL

if size(p,2) ~= 3
    error('input must be an Nx3 array');
end

rho = p(:,1);
phi = p(:,2);
z = p(:,3);

x = rho .* cos(phi);
y = rho .* sin(phi);

r = [x, y, z];
end