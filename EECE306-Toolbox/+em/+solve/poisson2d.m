function [V, info] = poisson2d(G, rho, epsr, mask, Vfix, opts)
%POISSON2D Solve div(eps grad V) = -rho on a uniform 2D grid.
%   [V,INFO] = EM.SOLVE.POISSON2D(G,RHO,EPSR,MASK,VFIX,OPTS) returns an
%   Ny-by-Nx potential. G comes from EM.SOLVE.GRID2D. RHO is a scalar,
%   Ny-by-Nx array, or function handle @(X,Y). EPSR is a positive scalar
%   or Ny-by-Nx array and is uniform for the Lab 7 stencil. MASK is true
%   where V is fixed to VFIX.
%
%   OPTS.method is 'jacobi', 'gs', 'sor', or 'direct'. OPTS.tol,
%   OPTS.maxit, and OPTS.omega default to 1e-8, 20000, and 1.9. Boundary
%   nodes outside MASK use homogeneous Neumann conditions. INFO contains
%   iter, resHist, time, and method.
%
%   See also EM.SOLVE.GRID2D, EM.TEST.MMS2D.

if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts,'method'), opts.method = 'sor'; end
if ~isfield(opts,'tol'),    opts.tol = 1e-8; end
if ~isfield(opts,'maxit'),  opts.maxit = 20000; end
if ~isfield(opts,'omega'),  opts.omega = 1.9; end
validateOptions(opts);

[Ny, Nx] = validateGrid(G);
hx2 = G.hx^2;
hy2 = G.hy^2;
e0 = em.const.eps0();

R = expand(rho, Ny, Nx, G, 'rho');
E = expand(epsr, Ny, Nx, G, 'epsr');

if any(E(:) <= 0)
    error('em:solve:poisson2d:InvalidPermittivity', ...
        'epsr values must be positive.');
end

if ~(islogical(mask) && isequal(size(mask), [Ny Nx]))
    error('em:solve:poisson2d:InvalidMask', ...
        'mask must be a logical Ny-by-Nx array.');
end

if ~(isnumeric(Vfix) && isreal(Vfix) && isequal(size(Vfix), [Ny Nx]) ...
        && all(isfinite(Vfix(:))))
    error('em:solve:poisson2d:InvalidFixedValues', ...
        'Vfix must be a finite real Ny-by-Nx array.');
end

% Right-hand side of the stencil equation, Laplacian(V) = S.
S = -R ./ (e0 * E);

% Five-point stencil for unequal x and y spacing.
cx = 1 / hx2;
cy = 1 / hy2;
cc = -2 * (cx + cy);

t0 = tic;
method = lower(opts.method);

switch method
    case 'direct'
        V = solveDirect(Ny, Nx, cx, cy, cc, S, mask, Vfix);
        info = struct('iter', 1, 'resHist', [], 'time', toc(t0), ...
            'method', 'direct');

    case {'jacobi', 'gs', 'sor'}
        [V, it, resHist] = solveIterative( ...
            Ny, Nx, cx, cy, cc, S, mask, Vfix, opts);

        info = struct('iter', it, 'resHist', resHist, 'time', toc(t0), ...
            'method', method);

    otherwise
        error('em:solve:poisson2d:BadMethod', ...
            'Unknown method %s.', opts.method);
end
end

function [V, it, resHist] = solveIterative( ...
        Ny, Nx, cx, cy, cc, S, mask, Vfix, opts)

V = zeros(Ny, Nx);
V(mask) = Vfix(mask);

jacobi = strcmpi(opts.method, 'jacobi');
sor = strcmpi(opts.method, 'sor');

w = opts.omega;
resHist = zeros(opts.maxit, 1);

for it = 1:opts.maxit

    Vold = V;

    for j = 1:Nx
        for i = 1:Ny

            if mask(i,j)
                continue
            end

            [iu, id, jl, jr] = neigh(i, j, Ny, Nx);

            if jacobi
                vu = Vold(iu,j);
                vd = Vold(id,j);
                vl = Vold(i,jl);
                vr = Vold(i,jr);
            else
                vu = V(iu,j);
                vd = V(id,j);
                vl = V(i,jl);
                vr = V(i,jr);
            end

            vnew = (S(i,j) - cx*(vl + vr) - cy*(vu + vd)) / cc;

            if sor
                V(i,j) = V(i,j) + w * (vnew - V(i,j));
            else
                V(i,j) = vnew;
            end
        end
    end

    r = norm(V(:) - Vold(:), inf);
    resHist(it) = r;

    if r < opts.tol
        break
    end
end

resHist = resHist(1:it);
end

function V = solveDirect(Ny, Nx, cx, cy, cc, S, mask, Vfix)

n = Ny * Nx;
idx = @(i,j) (j - 1) * Ny + i;

I = zeros(5*n,1);
J = I;
X = I;
b = zeros(n,1);
k = 0;

for j = 1:Nx
    for i = 1:Ny

        p = idx(i,j);

        if mask(i,j)
            k = k + 1;
            I(k) = p;
            J(k) = p;
            X(k) = 1;
            b(p) = Vfix(i,j);
            continue
        end

        [iu, id, jl, jr] = neigh(i, j, Ny, Nx);

        k = k + 1;
        I(k) = p;
        J(k) = p;
        X(k) = cc;

        k = k + 1;
        I(k) = p;
        J(k) = idx(i,jr);
        X(k) = cx;

        k = k + 1;
        I(k) = p;
        J(k) = idx(i,jl);
        X(k) = cx;

        k = k + 1;
        I(k) = p;
        J(k) = idx(iu,j);
        X(k) = cy;

        k = k + 1;
        I(k) = p;
        J(k) = idx(id,j);
        X(k) = cy;

        b(p) = S(i,j);
    end
end

A = sparse(I(1:k), J(1:k), X(1:k), n, n);

V = reshape(A \ b, Ny, Nx);
end

function [iu, id, jl, jr] = neigh(i, j, Ny, Nx)

if i == 1
    iu = 2;
    id = 2;
elseif i == Ny
    iu = Ny - 1;
    id = Ny - 1;
else
    iu = i + 1;
    id = i - 1;
end

if j == 1
    jr = 2;
    jl = 2;
elseif j == Nx
    jr = Nx - 1;
    jl = Nx - 1;
else
    jr = j + 1;
    jl = j - 1;
end

end

function M = expand(value, Ny, Nx, G, name)

if isa(value, 'function_handle')
    M = value(G.X, G.Y);
elseif isscalar(value)
    M = repmat(value, Ny, Nx);
else
    M = value;
end

if ~(isnumeric(M) && isreal(M) && isequal(size(M), [Ny Nx]) ...
        && all(isfinite(M(:))))
    error('em:solve:poisson2d:InvalidInput', ...
        '%s must produce a finite real Ny-by-Nx array.', name);
end
end

function [Ny, Nx] = validateGrid(G)

required = {'X','Y','hx','hy'};

if ~isstruct(G) || ~all(isfield(G, required))
    error('em:solve:poisson2d:InvalidGrid', ...
        'G must be a grid struct from em.solve.grid2d.');
end

[Ny, Nx] = size(G.X);

if ~(isnumeric(G.X) && isreal(G.X) && ...
        isnumeric(G.Y) && isreal(G.Y) && ...
        isequal(size(G.Y), [Ny Nx]) && ...
        Ny >= 2 && Nx >= 2 && ...
        all(isfinite(G.X(:))) && ...
        all(isfinite(G.Y(:))) && ...
        isnumeric(G.hx) && isscalar(G.hx) && ...
        isfinite(G.hx) && G.hx > 0 && ...
        isnumeric(G.hy) && isscalar(G.hy) && ...
        isfinite(G.hy) && G.hy > 0)

    error('em:solve:poisson2d:InvalidGrid', ...
        'G must contain finite Ny-by-Nx coordinates and positive spacings.');
end
end

function validateOptions(opts)

if ~(isstruct(opts) && isscalar(opts))
    error('em:solve:poisson2d:InvalidOptions', ...
        'opts must be a scalar struct.');
end

if ~(ischar(opts.method) && isrow(opts.method))
    error('em:solve:poisson2d:InvalidMethod', ...
        'opts.method must be text.');
end

if ~(isnumeric(opts.tol) && isreal(opts.tol) && ...
        isscalar(opts.tol) && isfinite(opts.tol) && opts.tol > 0)

    error('em:solve:poisson2d:InvalidTolerance', ...
        'opts.tol must be a finite positive scalar.');
end

if ~(isnumeric(opts.maxit) && isreal(opts.maxit) && ...
        isscalar(opts.maxit) && isfinite(opts.maxit) && ...
        opts.maxit == floor(opts.maxit) && opts.maxit >= 1)

    error('em:solve:poisson2d:InvalidIterationLimit', ...
        'opts.maxit must be a positive integer scalar.');
end

if ~(isnumeric(opts.omega) && isreal(opts.omega) && ...
        isscalar(opts.omega) && isfinite(opts.omega) && ...
        opts.omega > 0 && opts.omega < 2)

    error('em:solve:poisson2d:InvalidOmega', ...
        'opts.omega must be between 0 and 2.');
end
end
