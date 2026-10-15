% test_lab07.m  Two-dimensional Poisson solver and manufactured solutions.

e0 = em.const.eps0();

Vex = @(X,Y) sin(pi*X).*sin(pi*Y);
rho = em.test.mms2d(Vex, 1);

% Pointwise Laplacian check.
p = [0.2 0.3 0.4; -0.5 0.1 0.7];

fpoly = @(r) r(:,1).^2 + 2*r(:,2).^2 + 3*r(:,3).^2;

expectedL = 12 * ones(size(p,1),1);

em.test.assertClose( ...
    em.op.laplacian(fpoly,p,1e-3), ...
    expectedL, ...
    1e-8, ...
    'laplacian');

% 1. Manufactured solution error below 1e-3 on a 101 x 101 grid.
n = 101;

G = em.solve.grid2d([0 1],[0 1],n,n);

mk = false(n,n);
mk([1 end],:) = true;
mk(:,[1 end]) = true;

V = em.solve.poisson2d(G,rho,1,mk,zeros(n,n), ...
    struct('method','direct'));

err101 = max(abs(V(:)-Vex(G.X(:),G.Y(:))));

assert(err101 < 1e-3, ...
    '101-by-101 manufactured solution error');

% 2. Observed convergence order between 1.8 and 2.2.
Ns = [21 41 81 161]';
errs = zeros(size(Ns));

for k = 1:numel(Ns)

    n = Ns(k);

    Gn = em.solve.grid2d([0 1],[0 1],n,n);

    mk = false(n,n);
    mk([1 end],:) = true;
    mk(:,[1 end]) = true;

    Vn = em.solve.poisson2d(Gn,rho,1,mk,zeros(n,n), ...
        struct('method','direct'));

    errs(k) = max(abs(Vn(:)-Vex(Gn.X(:),Gn.Y(:))));
end

pOrder = em.test.convergence( ...
    @(N) errs(Ns==N), Ns, struct('plot',false));

assert(pOrder > 1.8 && pOrder < 2.2, ...
    'manufactured solution convergence order');

% 3 and 4. All methods agree and SOR < GS < Jacobi in iteration count.
n = 21;

G = em.solve.grid2d([0 1],[0 1],n,n);

mk = false(n,n);
mk([1 end],:) = true;
mk(:,[1 end]) = true;

Vfix = zeros(n,n);

Vd = em.solve.poisson2d(G,rho,1,mk,Vfix, ...
    struct('method','direct'));

methods = {'jacobi','gs','sor'};
iterations = zeros(1,3);

for k = 1:3

    opts = struct( ...
        'method',methods{k}, ...
        'tol',1e-9, ...
        'maxit',20000);

    if strcmp(methods{k},'sor')
        opts.omega = 2/(1+sin(pi/(n-1)));
    end

    [Vi,info] = em.solve.poisson2d(G,rho,1,mk,Vfix,opts);

    assert(max(abs(Vi(:)-Vd(:))) < 2e-7, ...
        'iterative and direct methods disagree');

    iterations(k) = info.iter;
end

assert(iterations(3) < iterations(2) && ...
       iterations(2) < iterations(1), ...
    'SOR must beat Gauss-Seidel, which must beat Jacobi');

% 5. Pure Laplace problem gives a linear potential profile.
n = 41;

G = em.solve.grid2d([0 1],[0 1],n,n);

mk = false(n,n);
mk(:,1) = true;
mk(:,end) = true;

Vfix = zeros(n,n);
Vfix(:,1) = 1;

VL = em.solve.poisson2d(G,0,1,mk,Vfix, ...
    struct('method','direct'));

ramp = repmat(linspace(1,0,n),n,1);

assert(max(abs(VL(:)-ramp(:))) < 1e-11, ...
    'linear Laplace profile');

% 6. No interior maximum or minimum in a source-free region.
inner = VL(2:end-1,2:end-1);

assert(min(inner(:)) > 0 && max(inner(:)) < 1, ...
    'interior maximum or minimum detected');

% 7. Gauss's law holds for the solver output.
n = 101;

G = em.solve.grid2d([-1 1],[-1 1],n,n);

mk = false(n,n);
mk([1 end],:) = true;
mk(:,[1 end]) = true;

rho0 = 1e-8;

rhoB = @(X,Y) rho0 * ...
    (abs(X)<0.25 & abs(Y)<0.25);

VB = em.solve.poisson2d(G,rhoB,1,mk,zeros(n,n), ...
    struct('method','direct'));

Vh = @(r) interp2(G.X,G.Y,VB, ...
    r(:,1),r(:,2),'linear');

Dh = @(r) -e0 * em.op.grad(Vh,r,G.hx/2);

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

    flux = flux + em.field.flux( ...
        Dh, faces{k}, ...
        [-a a], [-z z], ...
        20, 20, 'gauss');

end

Qenc = rho0 * 0.5 * 0.5 * (2*z);

assert(abs(flux-Qenc)/Qenc < 0.05, ...
    'Gauss law for solver output');

fprintf('  Lab 7: MMS error %.3e, order %.3f\n', ...
    err101, pOrder);

fprintf('  Jacobi / GS / SOR iterations %d / %d / %d\n', ...
    iterations);

fprintf('  Gauss flux %.4e C, enclosed %.4e C\n', ...
    flux, Qenc);

fprintf('  test_lab07 checks complete\n');
