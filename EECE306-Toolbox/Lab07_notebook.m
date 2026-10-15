%% EECE 306 Lab 7 notebook. Writing Your Own Poisson Solver
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore
%
% Run |publish('Lab07_notebook.m')| from the toolbox root, then print the
% HTML to PDF. A full run may take several minutes because the iterative
% solvers update one grid point at a time.

%% Setup
clear; close all;
e0 = em.const.eps0();

%% The manufactured solution
% The exact potential is chosen first. Its Laplacian gives the charge
% density, so the numerical solution can be compared with a known answer.
Vex = @(X,Y) sin(pi*X) .* sin(pi*Y);
rho = em.test.mms2d(Vex, 1);

n = 101;
G = em.solve.grid2d([0 1], [0 1], n, n);

mask = false(n,n);
mask(1,:) = true;
mask(end,:) = true;
mask(:,1) = true;
mask(:,end) = true;

Vfix = zeros(n,n);

[Vs, info] = em.solve.poisson2d(G, rho, 1, mask, Vfix, ...
    struct('method','direct'));

err101 = max(abs(Vs(:) - Vex(G.X(:), G.Y(:))));

fprintf('MMS max error on the 101 x 101 grid = %.3e\n', err101);

%% Convergence under grid refinement
Ns = [21 41 81 161]';
errs = zeros(size(Ns));

for j = 1:numel(Ns)

    n = Ns(j);

    Gn = em.solve.grid2d([0 1], [0 1], n, n);

    mk = false(n,n);
    mk(1,:) = true;
    mk(end,:) = true;
    mk(:,1) = true;
    mk(:,end) = true;

    Vn = em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), ...
        struct('method','direct'));

    errs(j) = max(abs(Vn(:) - Vex(Gn.X(:), Gn.Y(:))));
end

disp(table(Ns, errs, 'VariableNames', {'N','MaxError'}));

p = em.test.convergence(@(N) errs(Ns == N), Ns);

fprintf('observed order p = %.2f  (expect near 2 for the five point stencil)\n', p);

%% Three iterations and a direct solve
n = 41;

Gn = em.solve.grid2d([0 1], [0 1], n, n);

mk = false(n,n);
mk(1,:) = true;
mk(end,:) = true;
mk(:,1) = true;
mk(:,end) = true;

Vf = zeros(n,n);

methods = {'jacobi', 'gs', 'sor'};

figure;
hold on

for k = 1:3

    o = struct('method', methods{k}, ...
        'tol', 1e-6, 'maxit', 20000);

    if strcmp(methods{k}, 'sor')
        o.omega = 2/(1+sin(pi/(n-1)));
    end

    [Vk, ik] = em.solve.poisson2d(Gn, rho, 1, mk, Vf, o);

    semilogy(ik.resHist)

    fprintf('%-7s iterations = %5d\n', methods{k}, ik.iter);
end

grid on
legend(methods)
xlabel('iteration')
ylabel('residual')
title('Three iterative methods on the same 41 x 41 problem')

%% The omega scan on two grids
% We scan the same range on two grids to find where SOR needs the fewest
% iterations and to see whether the best value changes with grid size.

figure;
hold on

for n = [31 61]

    Gn = em.solve.grid2d([0 1], [0 1], n, n);

    mk = false(n,n);
    mk(1,:) = true;
    mk(end,:) = true;
    mk(:,1) = true;
    mk(:,end) = true;

    om = linspace(1.0, 1.99, 20);
    it = zeros(size(om));

    for j = 1:numel(om)

        [~, ij] = em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), ...
            struct('method','sor', ...
            'tol',1e-4, ...
            'maxit',20000, ...
            'omega',om(j)));

        it(j) = ij.iter;
    end

    plot(om, it, 'o-')

    [~, jb] = min(it);

    fprintf('n = %3d   best omega about %.2f\n', n, om(jb));
end

grid on
legend('31 x 31', '61 x 61')
xlabel('omega')
ylabel('iterations to tolerance')
title('SOR iteration count against omega')

%% Timing, iterative against direct
Nt = [41 61 81 101];
td = zeros(size(Nt));
ts = zeros(size(Nt));

for k = 1:numel(Nt)

    n = Nt(k);

    Gn = em.solve.grid2d([0 1], [0 1], n, n);

    mk = false(n,n);
    mk(1,:) = true;
    mk(end,:) = true;
    mk(:,1) = true;
    mk(:,end) = true;

    tic
    em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), ...
        struct('method','direct'));
    td(k) = toc;

    tic
    em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), ...
        struct('method','sor', ...
        'tol',1e-6, ...
        'maxit',40000, ...
        'omega',2/(1+sin(pi/(n-1)))));
    ts(k) = toc;

    fprintf('N = %3d   direct %.3f s   SOR %.3f s\n', ...
        n, td(k), ts(k));
end

figure;
loglog(Nt, td, 'o-', Nt, ts, 's-');
grid on

xlabel('grid size N')
ylabel('wall time (s)')
legend('direct', 'SOR', 'location', 'northwest')
title('Cost against grid size')

crossIdx = find(ts < td, 1);

if isempty(crossIdx)
    fprintf('No timing crossover was observed through N = %d.\n', Nt(end));
else
    fprintf('First observed timing crossover is near N = %d.\n', ...
        Nt(crossIdx));
end

%% A Laplace sanity pair
% One edge is held at 1 and the opposite edge at 0. The other two edges
% have zero normal derivative, so the potential should form a straight ramp.

n = 41;

Gn = em.solve.grid2d([0 1], [0 1], n, n);

mk = false(n,n);
mk(:,1) = true;
mk(:,end) = true;

Vf = zeros(n,n);
Vf(:,1) = 1;

VL = em.solve.poisson2d(Gn, 0, 1, mk, Vf, ...
    struct('method','direct'));

ramp = repmat(linspace(1, 0, n), n, 1);

fprintf('max deviation from the linear ramp = %.3e\n', ...
    max(abs(VL(:) - ramp(:))));

inner = VL(2:end-1, 2:end-1);

fprintf('interior range [%.4f, %.4f] inside boundary range [0, 1] = %d\n', ...
    min(inner(:)), max(inner(:)), ...
    min(inner(:)) >= 0 && max(inner(:)) <= 1);

%% Gauss's law closes the loop
% A block of charge sits inside a grounded box. The potential is solved
% first, then D is recovered from the gradient. The flux through a closed
% surface should equal the enclosed charge.

n = 101;

Gn = em.solve.grid2d([-1 1], [-1 1], n, n);

mk = false(n,n);
mk(1,:) = true;
mk(end,:) = true;
mk(:,1) = true;
mk(:,end) = true;

rho0 = 1e-8;

rhoB = @(X,Y) rho0 * (abs(X) < 0.25 & abs(Y) < 0.25);

VB = em.solve.poisson2d(Gn, rhoB, 1, mk, zeros(n,n), ...
    struct('method','direct'));

Vh = @(r) interp2(Gn.X, Gn.Y, VB, ...
    r(:,1), r(:,2), 'linear');

Dh = @(r) -e0 * em.op.grad(Vh, r, Gn.hx/2);

a = 0.6;
z = 0.5;

faces = { ...
    @(u,v)[u v z+0*u], ...
    @(u,v)[v u -z+0*u], ...
    @(u,v)[v a+0*u u], ...
    @(u,v)[u -a+0*u v], ...
    @(u,v)[a+0*u u v], ...
    @(u,v)[-a+0*u v u] ...
    };

flux = 0;

for k = 1:6
    flux = flux + em.field.flux(Dh, faces{k}, ...
        [-a a], [-z z], 20, 20, 'gauss');
end

Qenc = rho0 * 0.5 * 0.5 * (2*z);

fprintf('flux of D through the surface = %.4e C\n', flux);
fprintf('enclosed charge               = %.4e C\n', Qenc);
fprintf('relative difference           = %.2e\n', ...
    abs(flux-Qenc)/Qenc);

%% Interpretation
% The manufactured solution converged at about second order, which matches
% what we expect from the five-point stencil. If the measured order had
% been closer to one, the first things to check would be the boundary
% indexing and the Neumann boundary conditions. Gauss-Seidel needed fewer
% sweeps than Jacobi, and SOR was faster when omega was chosen near its best
% value. The best omega changed when the grid size changed, so one value
% does not work best for every grid. The linear potential and Gauss's law
% checks gave separate ways to verify the boundary conditions and the
% electric field recovered from the solution.

%% Problems encountered
% One problem we ran into was handling the zero-derivative boundary. At
% first, a missing neighbor outside the grid was replaced with the boundary
% node itself. The solution still looked smooth, so the mistake was not
% obvious. We corrected it by mirroring the interior neighbor across the
% boundary, which gives the centered zero-derivative condition required by
% the solver.

%% Full test suite
runTests
