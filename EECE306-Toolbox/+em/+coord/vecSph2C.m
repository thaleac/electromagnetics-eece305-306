function A = vecSph2C(As, r)
%VECSPH2C Converting spherical vector components to Cartesian components.
% A = EM.COORD.VECSPH2C(AS,R) converts Nx3 spherical vector
% components AS at Cartesian positions R.

if ~isnumeric(As)
    error('em:coord:vecC2Cyl:InvalidVectorType', ...
        'Input A must be a numeric Nx3 array.');
end

if ~isnumeric(r)
    error('em:coord:vecC2Cyl:InvalidPositionType', ...
        'Input r must be a numeric Nx3 Cartesian position array.');
end

if size(As,2) ~= 3 || size(r,2) ~= 3
    error('em:coord:vecC2Cyl:InvalidSize', ...
        'A and r must each have exactly 3 columns.');
end

if size(As,1) ~= size(r,1)
    error('em:coord:vecC2Cyl:RowMismatch', ...
        'A and r must contain the same number of rows.');
end

x = r(:,1);
y = r(:,2);
z = r(:,3);

rho = sqrt(x.^2 + y.^2);

theta = atan2(rho, z);
phi = atan2(y,x);

Ar = As(:,1);
Atheta = As(:,2);
Aphi = As(:,3);

Ax = Ar .* sin(theta) .* cos(phi) ...
    + Atheta .* cos(theta) .* cos(phi) ...
    - Aphi .* sin(phi);

Ay = Ar .* sin(theta) .* sin(phi) ...
    + Atheta .* cos(theta) .* sin(phi) ...
    + Aphi .* cos(phi);

Az = Ar .* cos(theta) ...
    - Atheta .* sin(theta);

A = [Ax Ay Az];
end
