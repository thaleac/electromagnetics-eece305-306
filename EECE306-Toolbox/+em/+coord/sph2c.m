function r = sph2c(s)
%SPH2C - From spherical to cartesian coordinates
%   R = EM.COORD.SPH2C(S) converts the Nx3 array of spherical points
%   S = [r theta phi] into Cartesian coordinates R = [x y z]
%
%   Example:
%       r = em.coord.sph2c([1 pi/2 0]) % returns [1 0 0]
%
%   See also EM.COORD.C2SPH

if size(s,2) ~= 3
    error('input must be an Nx3 array');
end

radius = s(:,1);
theta = s(:,2);
phi = s(:,3);

x = radius .* sin(theta) .* cos(phi);
y = radius .* sin(theta) .* sin(phi);
z = radius .* cos(theta);

r = [x, y, z];

end