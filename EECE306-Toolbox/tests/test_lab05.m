% test_lab05.m  Lab 5 flux, geometry, and divergence checks.
Q = 3e-9;
s = em.src.pointCharge(Q,[0 0 0]);
FD = @(r) em.field.D(s,r);
sph = @(th,ph) 0.5*[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
sphV = @(th,ph) 0.5*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];

% 1. D has the specified scaling, dimensions, and C/m^2 units.
pts = [1 0 0; 0 2 0; 0 0 -3];
Dv = em.field.D(s,pts);
em.test.assertClose(Dv,Q/(4*pi)*pts./sum(pts.^2,2).^(3/2),1e-12,'D closed form');
em.test.assertClose(em.field.D(s,pts,2.5),2.5*Dv,1e-12,'epsr scaling');
assert(isequal(size(em.field.D(s,zeros(0,3))),[0 3]),'D empty Nx3');

% 2. The required 1e-6 enclosed test uses resolved Gauss quadrature.
Phi = em.field.flux(FD,sph,[0 pi],[0 2*pi],32,64,'gauss');
em.test.assertClose(Phi,Q,1e-6,'enclosed sphere');
PhiFine = em.field.flux(FD,sphV,[0 pi],[0 2*pi],48,96,'gauss');
em.test.assertClose(PhiFine,Phi,1e-8,'Gauss refinement');

% Midpoint remains the documented default; its error must converge as N^-2.
Ns = [30 60 120]; errs = zeros(size(Ns));
for j = 1:numel(Ns)
    pm = em.field.flux(FD,sphV,[0 pi],[0 2*pi],Ns(j),2*Ns(j));
    errs(j) = abs(pm-Q)/Q;
end
assert(all(errs(1:2)./errs(2:3) > 3.9) && all(errs(1:2)./errs(2:3) < 4.1), ...
    'Midpoint sphere flux must converge at second order.');
em.test.assertClose(em.field.flux(FD,sphV,[0 pi],[0 2*pi],60,120,'midpoint'), ...
    Q*(pi/(120*sin(pi/120))),1e-8,'midpoint analytical discrete sum');
pm960 = em.field.flux(FD,sphV,[0 pi],[0 2*pi],960,1920);
em.test.assertClose(pm960,Q,1e-6,'resolved midpoint sphere');

% 3. External charge: absolute tolerance is 1e-6 of a hemisphere flux.
sout = em.src.pointCharge(Q,[5 0 0]);
Fout = @(r) em.field.D(sout,r);
pr = em.field.flux(Fout,sph,[0 pi],[-pi/2 pi/2],60,120);
pl = em.field.flux(Fout,sph,[0 pi],[pi/2 3*pi/2],60,120);
zeroTol = 8e-18; % ~1e-6 of the 7.445e-12 C hemisphere magnitude.
pout = em.field.flux(Fout,sph,[0 pi],[0 2*pi],32,64,'gauss');
poutFine = em.field.flux(Fout,sph,[0 pi],[0 2*pi],48,96,'gauss');
em.test.assertClose([pout poutFine],0,zeroTol,'outside sphere');
assert(abs(pr+pl) < 4e-16,'Midpoint hemisphere cancellation');

% 4. Three closed surfaces, using the starter's outward parameter orders.
s3 = em.src.pointCharge(Q,[0 0.05 0.1]); F3 = @(r) em.field.D(s3,r);
faces = { @(u,v)[u v 0.5*ones(size(u))], @(u,v)[v u -0.5*ones(size(u))], ...
    @(u,v)[v 0.5*ones(size(u)) u], @(u,v)[u -0.5*ones(size(u)) v], ...
    @(u,v)[0.5*ones(size(u)) u v], @(u,v)[-0.5*ones(size(u)) v u] };
normals = [0 0 1;0 0 -1;0 1 0;0 -1 0;1 0 0;-1 0 0];
pc = 0; pcFine = 0;
for k = 1:6
    pc = pc + em.field.flux(F3,faces{k},[-0.5 0.5],[-0.5 0.5],24,24,'gauss');
    pcFine = pcFine + em.field.flux(F3,faces{k},[-0.5 0.5],[-0.5 0.5],40,40,'gauss');
    % Geometry independently: constant fields integrate each signed area.
    for j = 1:3
        ej = zeros(1,3); ej(j) = 1;
        pa = em.field.flux(@(r) repmat(ej,size(r,1),1),faces{k}, ...
            [-0.5 0.5],[-0.5 0.5],4,4);
        em.test.assertClose(pa,normals(k,j),1e-9,'cube face signed area');
    end
end
cap = @(u,v) [u.*cos(v) -u.*sin(v) zeros(size(u))];
ps = em.field.flux(F3,sph,[0 pi],[0 2*pi],32,64,'gauss');
ph = em.field.flux(F3,sph,[0 pi/2],[0 2*pi],32,64,'gauss') ...
    + em.field.flux(F3,cap,[0 0.5],[0 2*pi],32,64,'gauss');
phFine = em.field.flux(F3,sph,[0 pi/2],[0 2*pi],48,96,'gauss') ...
    + em.field.flux(F3,cap,[0 0.5],[0 2*pi],48,96,'gauss');
em.test.assertClose([ps pc ph],Q,1e-4,'three surfaces enclosed charge');
em.test.assertClose([pcFine phFine],[pc ph],1e-7,'piecewise surface refinement');
em.test.assertClose([pc ph],ps,1e-4,'shape independence');

% Geometry before fields: sphere area, cap sign/area, and volume theorem.
[~,wa] = em.quad.surf(sph,[0 pi],[0 2*pi],32,64,'gauss');
em.test.assertClose(sum(wa),pi,1e-8,'radius 0.5 sphere area');
pz = em.field.flux(@(r) repmat([0 0 1],size(r,1),1),cap, ...
    [0 0.5],[0 2*pi],16,32,'gauss');
em.test.assertClose(pz,-pi*0.5^2,1e-8,'cap outward signed area');
em.test.assertClose(em.field.flux(@(r) r,sph,[0 pi],[0 2*pi],32,64,'gauss'), ...
    3*(4*pi*0.5^3/3),1e-8,'sphere geometry divergence theorem');

% 5. Three charges inside, two outside: only enclosed charges contribute.
charges = [2 -1 4 7 -3]*1e-9;
locations = [0 0 0;0.12 -0.04 0.1;-0.2 0.1 0.05;2 0 0;-1 2 0];
sm = em.src.pointCharge(charges(1),locations(1,:));
for k = 2:5
    sm = em.src.merge(sm,em.src.pointCharge(charges(k),locations(k,:)));
end
pmerge = em.field.flux(@(r) em.field.D(sm,r),sph,[0 pi],[0 2*pi],40,80,'gauss');
em.test.assertClose(pmerge,sum(charges(1:3)),1e-6,'merged enclosed charges');
fprintf('  merged source: flux %.8e C, enclosed %.8e C\n',pmerge,sum(charges(1:3)));

% 6. Independence of three off-center positions.
for x = [0 0.15 0.3]
    sm = em.src.pointCharge(Q,[x 0.1 0]);
    ppos = em.field.flux(@(r) em.field.D(sm,r),sph,[0 pi],[0 2*pi],40,80,'gauss');
    em.test.assertClose(ppos,Q,1e-6,'off-center position');
end
flip = @(ph,th) 0.5*[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
pflip = em.field.flux(FD,flip,[0 2*pi],[0 pi],64,32,'gauss');
em.test.assertClose(pflip,-Phi,1e-8,'orientation reversal');

% 7. Central derivative direction/scaling, default h, and O(h^2) error.
Fpoly = @(r) [r(:,1).^3+7*r(:,2), r(:,2).^3+11*r(:,3), r(:,3).^3+13*r(:,1)];
r = [0.2 -0.4 0.6; -0.5 0.3 0.7]; exact = 3*sum(r.^2,2);
e1 = em.op.div(Fpoly,r,0.02)-exact;
e2 = em.op.div(Fpoly,r,0.01)-exact;
em.test.assertClose(e1,3*0.02^2,1e-8,'central derivative denominator');
em.test.assertClose(e1./e2,4,1e-7,'second-order step-size scaling');
em.test.assertClose(em.op.div(Fpoly,r),exact+3e-10,1e-9,'default step');
assert(isequal(size(em.op.div(Fpoly,r)),[2 1]),'div must return Nx1');
assert(isequal(size(em.op.div(Fpoly,zeros(0,3))),[0 1]),'div empty Nx1');

% 8. Volume-density refinement at fixed h, including multiple interior points.
rhov = 1e-6; hdiv = 0.05;
ball = @(u,v,w) [w.*sin(u).*cos(v) w.*sin(u).*sin(v) w.*cos(u)];
rInside = [0.1 0.05 0.1;-0.1 -0.05 -0.1;-0.12 0.08 0.03];
dstudy = zeros(3,3);
for k = 1:3
    n = [20 40 80]; n = n(k);
    sb = em.src.volCharge(rhov,ball,[0 pi],[0 2*pi],[0 0.5],n,2*n,n);
    dstudy(k,:) = em.op.div(@(r) em.field.D(sb,r),rInside,hdiv)';
    em.test.assertClose(sum(sb.q.*sb.w),rhov*4*pi*0.5^3/3,5e-4,'ball charge geometry');
end
em.test.assertClose(dstudy(2:3,:),rhov,0.05,'interior density at 128000 and 1024000 elements');
assert(max(abs(dstudy(3,:)-dstudy(2,:)))/rhov < 0.03,'interior refinement stability');
em.test.assertClose(dstudy(:,1),dstudy(:,2),1e-8,'inversion symmetry');
Einside = em.field.D(sb,rInside);
em.test.assertClose(max(em.vec.mag(Einside-rhov*rInside/3)./em.vec.mag(rhov*rInside/3)), ...
    0,0.05,'uniform ball analytical interior field');
em.test.assertClose(em.field.D(sb,[0 0 0]),zeros(1,3),1e-18,'ball center symmetry');

% The continuum derivative is rho, but the discrete point-element field has
% zero true divergence between elements: shrinking h need not improve it.
dtiny = em.op.div(@(r) em.field.D(sb,r),rInside(1,:),1e-5);
assert(abs(dtiny) < 0.01*rhov,'tiny h resolves charge-free element gaps');
fprintf('  density errors at 128000 / 1024000 elements: %.3f%% / %.3f%% max\n', ...
    100*max(abs(dstudy(2,:)-rhov))/rhov,100*max(abs(dstudy(3,:)-rhov))/rhov);

% 9. Charge-free divergence: finite-step error shrinks quadratically.
dout1 = em.op.div(FD,[2 1 1],0.05);
dout2 = em.op.div(FD,[2 1 1],0.025);
em.test.assertClose(dout1,0,2e-14,'charge-free divergence finite h');
em.test.assertClose(em.op.div(FD,[2 1 1]),0,1e-19,'charge-free divergence default h');
assert(abs(dout2) < 0.3*abs(dout1),'charge-free divergence improves with h halved');

% 10. Planar flux supports every rule and scalar-only/vectorized geometry.
plane = @(u,v) [u v 0];
for rule = {'midpoint','trapz','simpson','gauss'}
    pa = em.field.flux(@(r) repmat([0 0 2],size(r,1),1),plane,[0 2],[0 3],5,7,rule{1});
    em.test.assertClose(pa,12,1e-9,'planar flux field times area');
end

% 11. Invalid shapes and nonfinite fields must not silently pass.
badCalls = { @() em.field.D(s,[1 2]), @() em.field.D(s,[1 0 0],[1 2]), ...
    @() em.field.flux(FD,plane,[0 1 2],[0 1],4,4), ...
    @() em.field.flux(FD,plane,[0 1],[0 1],2.5,4), ...
    @() em.field.flux(@(r) ones(size(r,1),1),plane,[0 1],[0 1],4,4), ...
    @() em.op.div(FD,[1 2]), @() em.op.div(FD,[1 0 0],0), ...
    @() em.op.div(@(r) NaN(size(r)),[1 0 0]), ...
    @() em.field.flux(FD,@(u,v)[0 0 0],[0 1],[0 1],1,1) };
for k = 1:numel(badCalls)
    caught = false;
    try
        badCalls{k}();
    catch err
        caught = true;
        if k == numel(badCalls)
            assert(~isempty(strfind(err.message,'element 1')),'flux must name singular source element');
        end
    end
    assert(caught,'Invalid input or singularity must raise an error.');
end
pts10k = repmat([2 1 1],10000,1);
tic; d10k = em.op.div(FD,pts10k); elapsed = toc;
assert(isequal(size(d10k),[10000 1]) && elapsed < 2,'vectorized 10000-point divergence');
fprintf('  sphere / cube / closed hemisphere: %.8e / %.8e / %.8e C\n',ps,pc,ph);
fprintf('  sphere Gauss rel err %.2e; midpoint N=960 rel err %.2e\n',abs(Phi-Q)/Q,abs(pm960-Q)/Q);
disp('  test_lab05 checks complete');
