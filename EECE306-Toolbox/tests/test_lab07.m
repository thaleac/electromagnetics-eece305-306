% test_lab07.m  Two-dimensional Poisson solver and manufactured solutions.
e0 = em.const.eps0();
Vex = @(X,Y) sin(pi*X).*sin(pi*Y);
rho = em.test.mms2d(Vex, 1);

% Grid contract and independent pointwise Laplacian.
G = em.solve.grid2d([0 1], [0 2], 11, 21);
assert(isequal(size(G.X),[21 11]) && isequal(size(G.Y),[21 11]));
em.test.assertClose([G.hx G.hy], [0.1 0.1], 1e-14, 'grid spacing');
p = [0.2 0.3 0.4; -0.5 0.1 0.7];
fpoly = @(r) r(:,1).^2 + 2*r(:,2).^2 + 3*r(:,3).^2;
em.test.assertClose(em.op.laplacian(fpoly,p,1e-3),12,1e-8,'laplacian');

% Manufactured solution error and second-order convergence.
Ns = [21 41 81 101]'; errs = zeros(size(Ns));
for k = 1:numel(Ns)
    n = Ns(k); Gn = em.solve.grid2d([0 1],[0 1],n,n);
    mk = false(n,n); mk([1 end],:) = true; mk(:,[1 end]) = true;
    Vn = em.solve.poisson2d(Gn,rho,1,mk,0,struct('method','direct'));
    errs(k) = max(abs(Vn(:)-Vex(Gn.X(:),Gn.Y(:))));
end
assert(errs(end)<1e-3,'101-by-101 MMS error');
pOrder = em.test.convergence(@(N) errs(Ns==N),Ns,struct('plot',false));
assert(pOrder>1.8 && pOrder<2.2,'MMS convergence order');

% All four methods must agree, with the expected iteration ordering.
n = 21; Gi = em.solve.grid2d([0 1],[0 1],n,n);
mk = false(n,n); mk([1 end],:) = true; mk(:,[1 end]) = true;
[Vd,id] = em.solve.poisson2d(Gi,rho,1,mk,0,struct('method','direct'));
names = {'jacobi','gs','sor'}; iterations = zeros(1,3);
for k = 1:3
    opts = struct('method',names{k},'tol',1e-9,'maxit',20000);
    if strcmp(names{k},'sor'), opts.omega = 2/(1+sin(pi/(n-1))); end
    [Vi,ii] = em.solve.poisson2d(Gi,rho,1,mk,0,opts);
    assert(max(abs(Vi(:)-Vd(:)))<2e-7,'iterative/direct agreement');
    assert(ii.resHist(end)<opts.tol && numel(ii.resHist)==ii.iter);
    iterations(k) = ii.iter;
end
assert(iterations(3)<iterations(2) && iterations(2)<iterations(1), ...
    'SOR must beat Gauss-Seidel, which must beat Jacobi');
assert(id.iter==1 && isempty(id.resHist));

% Linear Laplace profile and the discrete maximum principle.
n = 41; GL = em.solve.grid2d([0 1],[0 1],n,n);
mk = false(n,n); mk(:,1) = true; mk(:,end) = true;
Vfix = zeros(n,n); Vfix(:,1) = 1;
VL = em.solve.poisson2d(GL,0,1,mk,Vfix,struct('method','direct'));
ramp = repmat(linspace(1,0,n),n,1);
assert(max(abs(VL(:)-ramp(:)))<1e-11,'linear Laplace profile');
inner = VL(2:end-1,2:end-1);
assert(min(inner(:))>0 && max(inner(:))<1,'discrete maximum principle');

% Gauss law for a field recovered from the solved potential.
n = 101; GB = em.solve.grid2d([-1 1],[-1 1],n,n);
mk = false(n,n); mk([1 end],:) = true; mk(:,[1 end]) = true;
rho0 = 1e-8; rhoB = @(X,Y) rho0*(abs(X)<0.25 & abs(Y)<0.25);
VB = em.solve.poisson2d(GB,rhoB,1,mk,0,struct('method','direct'));
Vh = @(r) interp2(GB.X,GB.Y,VB,r(:,1),r(:,2),'linear');
Dh = @(r) -e0*em.op.grad(Vh,r,GB.hx/2);
a = 0.6; z = 0.5;
faces = { @(u,v)[u v z+0*u], @(u,v)[v u -z+0*u], ...
    @(u,v)[v a+0*u u], @(u,v)[u -a+0*u v], ...
    @(u,v)[a+0*u u v], @(u,v)[-a+0*u v u] };
flux = 0;
for k = 1:6
    flux = flux + em.field.flux(Dh,faces{k},[-a a],[-z z],20,20,'gauss');
end
Qenc = rho0 * 0.5 * 0.5 * (2*z);
assert(abs(flux-Qenc)/Qenc<0.05,'Gauss law for solver output');

% Error contracts for the four new functions.
bad = {@() em.solve.grid2d([1 0],[0 1],3,3), ...
    @() em.op.laplacian(@(r) r(:,1),[1 2],1e-3), ...
    @() em.test.mms2d(Vex,0), ...
    @() em.solve.poisson2d(GL,0,1,false(size(GL.X)),0,struct('method','bad'))};
for k = 1:numel(bad)
    caught = false; try, bad{k}(); catch, caught = true; end; assert(caught);
end
fprintf('  Lab 7: MMS error %.3e, order %.3f\n',errs(end),pOrder);
fprintf('  Jacobi / GS / SOR iterations %d / %d / %d\n',iterations);
fprintf('  Gauss flux %.4e C, enclosed %.4e C\n',flux,Qenc);
fprintf('  test_lab07 checks complete\n');
