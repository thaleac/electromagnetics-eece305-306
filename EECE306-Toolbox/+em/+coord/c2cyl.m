function p = c2cyl(r)
%C2CYL - From cartesian to cylindrical coordinates
%   P = EM.COORD.C2CYL(R) converts the Nx3 array of Cartesian points
%   R = [x y z] into cylindrical coordinates [rho phi z]
%
%   Example:
%       p = em.coord.c2cyl([1 1 0]) % returns [sqrt(2)  pi/4  0]
%
%   See also EM.COORD.CYL2C

if size(r,2) ~= 3
    error('input must be an Nx3 array');
end
x = r(:,1);
y = r(:,2);
z = r(:,3);

rho = sqrt(x.^2 + y.^2);
phi = mod(atan2(y,x), 2*pi);
p = [rho, phi, z];
end