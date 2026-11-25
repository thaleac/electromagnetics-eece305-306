function s = c2sph(r)
%C2SPH - From cartesian to spherical coordinates
%   S = EM.COORD.C2SPH(R) converts the Nx3 array of Cartesian points
%   R = [x y z] into spherical coordinates S = [r theta phi]
%
%   Example:
%       s = em.coord.c2sph([1 1 1]) % returns [sqrt(3) atan2(sqrt(2),1) pi/4]
%
%   See also EM.COORD.SPH2C

if size(r,2) ~= 3
    error('input must be an Nx3 array');
end

x = r(:,1);
y = r(:,2);
z = r(:,3);

rho = sqrt(x.^2 + y.^2);
radius = sqrt(x.^2 + y.^2 + z.^2);
theta = atan2(rho, z);
phi = mod(atan2(y,x), 2*pi);

s = [radius, theta, phi];

end