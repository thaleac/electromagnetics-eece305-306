% test_lab04.m  Lab 4 quadrature and continuous-charge checks.
e0 = em.const.eps0();

% 1. Geometry weights recover length, area, and volume.
circ = @(t) [cos(t) sin(t) zeros(size(t))];
[~,wc] = em.quad.line(circ,[0 2*pi],400);
em.test.assertClose(sum(wc),2*pi,1e-6,'unit circle length');
sph = @(th,ph) [sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
[~,ws] = em.quad.surf(sph,[0 pi],[0 2*pi],32,32,'gauss');
em.test.assertClose(sum(ws),4*pi,1e-8,'unit sphere area');
ball = @(u,v,w) [w*sin(u)*cos(v) w*sin(u)*sin(v) w*cos(u)];
[~,wv] = em.quad.vol(ball,[0 pi],[0 2*pi],[0 1],16,16,8,'gauss');
em.test.assertClose(sum(wv),4*pi/3,1e-8,'unit ball volume');
sq = @(u,v) [u v 0];
[~,wq] = em.quad.surf(sq,[0 1],[0 1],20,20);
em.test.assertClose(sum(wq),1,1e-10,'unit square area');

% 2. Simpson and two-point Gauss integrate a cubic exactly.
exactCubic = 1/4 + 2/3 + 3/2 + 4;
[t,w] = em.quad.nodes(0,1,5,'simpson');
em.test.assertClose(sum((t.^3+2*t.^2+3*t+4).*w),exactCubic,1e-12, ...
    'Simpson cubic exactness');
[t,w] = em.quad.nodes(0,1,2,'gauss');
em.test.assertClose(sum((t.^3+2*t.^2+3*t+4).*w),exactCubic,1e-12, ...
    'Gauss cubic exactness');

% 3. Measured orders on the smooth integral of exp(t) from zero to one.
Iex = exp(1)-1;
Ns = [17 33 65 129]';
errMid = zeros(size(Ns));
errTrap = zeros(size(Ns));
errSimp = zeros(size(Ns));
for k = 1:numel(Ns)
    [t,w] = em.quad.nodes(0,1,Ns(k),'midpoint');
    errMid(k) = abs(sum(exp(t).*w)-Iex);
    [t,w] = em.quad.nodes(0,1,Ns(k),'trapz');
    errTrap(k) = abs(sum(exp(t).*w)-Iex);
    [t,w] = em.quad.nodes(0,1,Ns(k),'simpson');
    errSimp(k) = abs(sum(exp(t).*w)-Iex);
end
pMid = em.test.convergence(@(N) errMid(Ns==N),Ns,struct('plot',false));
pTrap = em.test.convergence(@(N) errTrap(Ns==N),Ns,struct('plot',false));
pSimp = em.test.convergence(@(N) errSimp(Ns==N),Ns,struct('plot',false));
assert(abs(pMid-2) < 0.1,'Midpoint observed order must be near 2.');
assert(abs(pTrap-2) < 0.1,'Trapezoid observed order must be near 2.');
assert(abs(pSimp-4) < 0.2,'Simpson observed order must be near 4.');

% 4. A finite line looks like its total point charge from far away.
lineC = @(t) [zeros(size(t)) zeros(size(t)) t];
sl = em.src.lineCharge(1e-9,lineC,[-0.5 0.5],800);
Eline = em.field.E(sl,[20 0 0]);
Epoint = em.field.E(em.src.pointCharge(1e-9,[0 0 0]),[20 0 0]);
assert(abs(Eline(1)-Epoint(1))/abs(Epoint(1)) < 0.01, ...
    'Far line-charge field must agree with the equivalent point charge.');

% 5. The sampled finite-line threshold is bracketed by adjacent L/d values.
Lbelow = 13.8950; Labove = 19.3070; d = 1;
sBelow = em.src.lineCharge(1e-9,lineC,[-Lbelow/2 Lbelow/2],800);
sAbove = em.src.lineCharge(1e-9,lineC,[-Labove/2 Labove/2],800);
rBelow = em.field.E(sBelow,[d 0 0]);
rAbove = em.field.E(sAbove,[d 0 0]);
rBelow = rBelow(1)/(1e-9/(2*pi*e0*d));
rAbove = rAbove(1)/(1e-9/(2*pi*e0*d));
assert(rBelow < 0.99 && rAbove > 0.99, ...
    'Finite-line 0.99 threshold must be bracketed by the reported samples.');

% 6. The sampled disk threshold exceeds 0.99 at a/h about 59.8.
hgt = 1; a = 59.8306;
disk = @(u,v) [u.*cos(v) u.*sin(v) zeros(size(u))];
sd = em.src.surfCharge(1e-9,disk,[0 a],[0 2*pi],80,80);
Ed = em.field.E(sd,[0 0 hgt]);
ratioD = Ed(3)/(1e-9/(2*e0));
assert(ratioD > 0.99,'Disk ratio must exceed 0.99 at reported a/h.');

% 7. An odd density integrates to zero on a symmetric line.
sn = em.src.lineCharge(@(r) r(:,3),lineC,[-1 1],400);
em.test.assertClose(sum(sn.q.*sn.w),0,1e-12, ...
    'odd line density net charge');
disp('  test_lab04 checks complete');
