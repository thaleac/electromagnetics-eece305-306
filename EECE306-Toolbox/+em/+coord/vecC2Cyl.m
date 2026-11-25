function Ac = vecC2Cyl(A, r)
%VECC2CYL Convert Cartesian vector components to cylindrical components.
% AC = EM.COORD.VECC2CYL(A,R) converts Nx3 Cartesisan vector
% components A at Cartesian positions R

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

phi = atan2(y,x);

Ax = A(:,1);
Ay = A(:,2);
Az = A(:,3);

Arho = Ax .* cos(phi) + Ay .* sin(phi);
Aphi = -Ax .* sin(phi) + Ay .* cos(phi);

Ac = [Arho Aphi Az];
end