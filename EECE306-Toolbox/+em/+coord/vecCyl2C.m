function A = vecCyl2C(Ac, r)
%VECCYL2C converting cylindrical vector components to Cartesian components.
%A = Em.COORD.VECCYL2C(AC,R) converts Nx3 cylindrical vector
% components AC at Cartesian postions R.

if ~isnumeric(Ac)
    error('em:coord:vecC2Sph:InvalidVectorType', ...
        'Input A must be a numeric Nx3 array.');
end

if ~isnumeric(r)
    error('em:coord:vecC2Sph:InvalidPositionType', ...
        'Input r must be a numeric Nx3 Cartesian position array.');
end

if size(Ac,2) ~= 3 || size(r,2) ~= 3
    error('em:coord:vecC2Sph:InvalidSize', ...
        'A and r must each have exactly 3 columns.');
end

if size(Ac,1) ~= size(r,1)
    error('em:coord:vecC2Sph:RowMismatch', ...
        'A and r must contain the same number of rows.');
end

x = r(:,1);
y = r(:,2);

phi = atan2(y,x);

Arho = Ac(:,1);
Aphi = Ac(:,2);
Az = Ac(:,3);

Ax = Arho .* cos(phi) - Aphi .* sin(phi);
Ay = Arho .* sin(phi) + Aphi .* cos(phi);

A = [Ax Ay Az];
end
