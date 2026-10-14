function G = grid2d(xlim, ylim, Nx, Ny)
%GRID2D Build a uniform two-dimensional finite-difference grid.
%   G = EM.SOLVE.GRID2D(XLIM,YLIM,NX,NY) returns row coordinate vectors
%   G.x and G.y, Ny-by-Nx meshgrid arrays G.X and G.Y, spacings G.hx and
%   G.hy, and the node counts G.Nx and G.Ny.
%
%   Example:
%       G = em.solve.grid2d([0 1], [0 2], 51, 81);
%
%   See also MESHGRID, EM.SOLVE.POISSON2D.

validateLimit(xlim, 'xlim');
validateLimit(ylim, 'ylim');
validateCount(Nx, 'Nx');
validateCount(Ny, 'Ny');

x = linspace(xlim(1), xlim(2), Nx);
y = linspace(ylim(1), ylim(2), Ny);
[X, Y] = meshgrid(x, y);
G = struct('x', x, 'y', y, 'X', X, 'Y', Y, ...
    'hx', x(2)-x(1), 'hy', y(2)-y(1), 'Nx', Nx, 'Ny', Ny);
end

function validateLimit(value, name)
if ~(isnumeric(value) && isreal(value) && isequal(size(value), [1 2]) ...
        && all(isfinite(value)) && value(2) > value(1))
    error('em:solve:grid2d:InvalidLimit', ...
        '%s must be a finite increasing 1x2 array.', name);
end
end

function validateCount(value, name)
if ~(isnumeric(value) && isreal(value) && isscalar(value) ...
        && isfinite(value) && value == floor(value) && value >= 2)
    error('em:solve:grid2d:InvalidCount', ...
        '%s must be an integer scalar of at least 2.', name);
end
end
