function h = quiver2(F, xBounds, yBounds, N, opts)
%QUIVER2 Plot a two-dimensional sample of a vector field.
%   H = EM.VIZ.QUIVER2(F, XLIM, YLIM, N, OPTS) samples field function F
%   on an NxN grid in the z=0 plane and draws its x-y components.
%   OPTS supports normalize, logscale, and title.
%
%   Example:
%       s = em.src.pointCharge(1e-9, [0 0 0]);
%       em.viz.quiver2(@(r) em.field.E(s,r), [-1 1], [-1 1], 16, ...
%           struct('normalize', true, 'title', 'Point-charge field'));
%
%   See also QUIVER, EM.FIELD.E.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isa(F, 'function_handle')
    error('em:viz:quiver2:InvalidField', ...
        'F must be a field function handle.');
end
if ~(isnumeric(xBounds) && isreal(xBounds) && isequal(size(xBounds), [1 2]) ...
        && all(isfinite(xBounds)) && xBounds(1) < xBounds(2))
    error('em:viz:quiver2:InvalidXLim', ...
        'xlim must be a finite increasing 1x2 array.');
end
if ~(isnumeric(yBounds) && isreal(yBounds) && isequal(size(yBounds), [1 2]) ...
        && all(isfinite(yBounds)) && yBounds(1) < yBounds(2))
    error('em:viz:quiver2:InvalidYLim', ...
        'ylim must be a finite increasing 1x2 array.');
end
if ~(isnumeric(N) && isscalar(N) && isfinite(N) && N >= 2 && N == floor(N))
    error('em:viz:quiver2:InvalidN', ...
        'N must be an integer scalar greater than or equal to 2.');
end
if ~isstruct(opts)
    error('em:viz:quiver2:InvalidOptions', 'opts must be a struct.');
end

normalize = getOption(opts, 'normalize', false);
logscale = getOption(opts, 'logscale', false);
plotTitle = getOption(opts, 'title', 'Vector field');
if ~(islogical(normalize) && isscalar(normalize))
    error('em:viz:quiver2:InvalidNormalize', ...
        'opts.normalize must be a logical scalar.');
end
if ~(islogical(logscale) && isscalar(logscale))
    error('em:viz:quiver2:InvalidLogscale', ...
        'opts.logscale must be a logical scalar.');
end
if ~(ischar(plotTitle) || (isstring(plotTitle) && isscalar(plotTitle)))
    error('em:viz:quiver2:InvalidTitle', ...
        'opts.title must be character text.');
end

[X, Y] = meshgrid(linspace(xBounds(1), xBounds(2), N), ...
                  linspace(yBounds(1), yBounds(2), N));
points = [X(:) Y(:) zeros(numel(X), 1)];
Ev = F(points);
if ~isnumeric(Ev) || ~isequal(size(Ev), [numel(X) 3])
    error('em:viz:quiver2:InvalidFieldOutput', ...
        'F must return an Nx3 numeric field for an Nx3 input.');
end

mag = sqrt(sum(Ev.^2, 2));
if normalize
    nz = mag > 0;
    Ev(nz, :) = Ev(nz, :) ./ mag(nz);
elseif logscale
    nz = mag > 0;
    Ev(nz, :) = Ev(nz, :) .* (log1p(mag(nz)) ./ mag(nz));
end

U = reshape(Ev(:,1), size(X));
V = reshape(Ev(:,2), size(Y));
h = quiver(X, Y, U, V);
axis equal;
xlim(xBounds);
ylim(yBounds);
xlabel('x (m)');
ylabel('y (m)');
title(plotTitle);
grid on;
end

function value = getOption(opts, name, defaultValue)
if isfield(opts, name)
    value = opts.(name);
else
    value = defaultValue;
end
end
