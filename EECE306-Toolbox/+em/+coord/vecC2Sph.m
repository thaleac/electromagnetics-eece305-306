function As = vecC2Sph(A, r)
%VECC2SPH Converting Cartesian vector components to spherical components.
% AS = EM.COORD.VECC2SPH(A,R) converts Nx3 Cartesian vector
% components A at Cartesian positions R.

if ~isnumeric(A)
    error('em:coord:vecC2Sph:InvalidVectorType', ...
        'Input A must be a numeric Nx3 array.');
end

if ~isnumeric(r)
    error('em:coord:vecC2Sph:InvalidPositionType', ...
        'Input r must be a numeric Nx3 Cartesian position array.');
end

if size(A,2) ~= 3 || size(r,2) ~= 3
    error('em:coord:vecC2Sph:InvalidSize', ...
        'A and r must each have exactly 3 columns.');
end

if size(A,1) ~= size(r,1)
    error('em:coord:vecC2Sph:RowMismatch', ...
        'A and r must contain the same number of rows.');
end

x = r(:,1);
y = r(:,2);
z = r(:,3);

rho = sqrt(x.^2 + y.^2);

theta = atan2(rho, z);
phi = atan2(y,x);

Ax = A(:,1);
Ay = A(:,2);
Az = A(:,3);

Ar = Ax .* sin(theta) .* cos(phi) ...
   + Ay .* sin(theta) .* sin(phi) ...
   + Az .* cos(theta);


Atheta = Ax .* cos(theta) .* cos(phi) ...
       + Ay .* cos(theta) .* sin(phi)...
       - Az .* sin(theta);

Aphi = -Ax .* sin(phi)...
     + Ay .* cos(phi);

As = [Ar Atheta Aphi];
end

