%% EECE 306 Lab 7 notebook. Writing Your Own Poisson Solver
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore

%% Setup
clear; close all;
e0 = em.const.eps0();

Vex = @(X,Y) sin(pi*X).*sin(pi*Y);
rho = em.test.mms2d(Vex, 1);

%% Manufactured solution and convergence
n = 101;
G = em.solve.grid2d([0 1], [0 1], n, n);
mask = false(n,n);
mask([1 end],:) = true;
mask(:,[1 end]) = true;

V = em.solve.poisson2d(G, rho, 1, mask, zeros(n,n), ...
    struct('method','direct'));

err101 = max(abs(V(:) - Vex(G.X(:),G.Y(:))));
fprintf('MMS max error on the 101 x 101 grid = %.3e\n', err101);

Ns = [21 41 81 161]';
errs = zeros(size(Ns));

for j = 1:numel(Ns)
    n = Ns(j);
    Gn = em.solve.grid2d([0 1], [0 1], n, n);
    mk = false(n,n);
    mk([1 end],:) = true;
    mk(:,[1 end]) = true;

    Vn = em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), ...
        struct('method','direct'));

    errs(j) = max(abs(Vn(:) - Vex(Gn.X(:),Gn.Y(:))));
end

disp(table(Ns, errs, 'VariableNames', {'N','MaxError'}));

p = em.test.convergence(@(N) errs(Ns == N), Ns);
fprintf('observed order p = %.2f\n', p);

%% Three iterative methods
n = 41;
G = em.solve.grid2d([0 1], [0 1], n, n);
mask = false(n,n);
mask([1 end],:) = true;
mask(:,[1 end]) = true;

methods = {'jacobi','gs','sor'};

figure;
hold on

for k = 1:3
    opts = struct('method',methods{k}, ...
        'tol',1e-6,'maxit',20000);

    if strcmp(methods{k},'sor')
        opts.omega = 2/(1+sin(pi/(n-1)));
    end

    [~,info] = em.solve.poisson2d(G, rho, 1, mask, ...
        zeros(n,n), opts);

    semilogy(info.resHist)
    fprintf('%-7s iterations = %d\n', ...
        methods{k}, info.iter);
end

grid on
legend(methods)
xlabel('iteration')
ylabel('residual')
title('Three iterative methods on the same 41 x 41 problem')

%% Omega scan on two grids
figure;
hold on

for n = [31 61]
    G = em.solve.grid2d([0 1], [0 1], n, n);
    mask = false(n,n);
    mask([1 end],:) = true;
    mask(:,[1 end]) = true;

    omega = linspace(1.0,1.99,20);
    iterations = zeros(size(omega));

    for j = 1:numel(omega)
        [~,info] = em.solve.poisson2d(G, rho, 1, mask, ...
            zeros(n,n), ...
            struct('method','sor','tol',1e-4, ...
            'maxit',20000,'omega',omega(j)));

        iterations(j) = info.iter;
    end

    plot(omega,iterations,'o-')

    [~,best] = min(iterations);
    fprintf('n = %d   best omega about %.2f\n', ...
        n, omega(best));
end

grid on
legend('31 x 31','61 x 61')
xlabel('omega')
ylabel('iterations to tolerance')
title('SOR iteration count against omega')

%% Timing, iterative against direct
Nt = [41 61 81 101];
td = zeros(size(Nt));
ts = zeros(size(Nt));

for k = 1:numel(Nt)
    n = Nt(k);

    G = em.solve.grid2d([0 1], [0 1], n, n);
    mask = false(n,n);
    mask([1 end],:) = true;
    mask(:,[1 end]) = true;

    tic
    em.solve.poisson2d(G, rho, 1, mask, zeros(n,n), ...
        struct('method','direct'));
    td(k) = toc;

    tic
    em.solve.poisson2d(G, rho, 1, mask, zeros(n,n), ...
        struct('method','sor','tol',1e-6, ...
        'maxit',40000, ...
        'omega',2/(1+sin(pi/(n-1)))));
    ts(k) = toc;

    fprintf('N = %d   direct %.3f s   SOR %.3f s\n', ...
        n, td(k), ts(k));
end

figure;
loglog(Nt,td,'o-',Nt,ts,'s-');
grid on
xlabel('grid size N')
ylabel('wall time (s)')
legend('direct','SOR','location','northwest')
title('Cost against grid size')

crossIdx = find(ts < td,1);

if isempty(crossIdx)
    fprintf('No timing crossover was observed through N = %d.\n', ...
        Nt(end));
else
    fprintf('First observed timing crossover is near N = %d.\n', ...
        Nt(crossIdx));
end

%% Interpretation
% The manufactured solution converged at about second order, which matches
% the five-point stencil. If the measured order had been closer to one, the
% first things to check would be the boundary indexing and the Neumann
% boundary conditions. Gauss-Seidel needed fewer sweeps than Jacobi, while
% SOR was faster when omega was chosen near its best value. The best omega
% changed as the grid became finer, showing that one value does not work
% best for every grid.

%% Problems encountered
% One problem we ran into was handling the zero-derivative boundary. At
% first, a missing neighbor outside the grid was replaced with the boundary
% node itself. The solution still looked smooth, so the mistake was not
% obvious. We corrected it by mirroring the interior neighbor across the
% boundary, which gives the centered zero-derivative condition required by
% the solver.

%% Full test suite
runTests
close all