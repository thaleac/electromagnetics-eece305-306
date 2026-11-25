% test_lab03.m  Lab 3 point-charge, field, and visualization checks.
rng(306);
e0 = em.const.eps0();

% 1. Closed form at ten radii.
Q = 2e-9;
s = em.src.pointCharge(Q, [0 0 0]);
radii = (1:10)';
Ev = em.field.E(s, [radii zeros(10,2)]);
Eref = Q ./ (4*pi*e0*radii.^2);
em.test.assertClose(Ev(:,1), Eref, 1e-12, 'point-charge closed form');
em.test.assertClose(Ev(:,2:3), zeros(10,2), 1e-12, 'point-charge transverse field');

% 2. Positive fields point outward and negative fields point inward.
dirs = em.vec.unit([1 2 3; -2 1 4; 3 -4 2]);
sp = em.src.pointCharge(1e-9, [0 0 0]);
sn = em.src.pointCharge(-1e-9, [0 0 0]);
assert(all(sum(em.field.E(sp, dirs).*dirs,2) > 0), ...
    'Positive point-charge field must point radially outward.');
assert(all(sum(em.field.E(sn, dirs).*dirs,2) < 0), ...
    'Negative point-charge field must point radially inward.');

% 3. Symmetry on the perpendicular bisector plane.
q1 = em.src.pointCharge(1e-9, [-0.1 0 0]);
q2 = em.src.pointCharge(1e-9, [ 0.1 0 0]);
s2 = em.src.merge(q1, q2);
pts = [zeros(200,1) randn(200,2)];
Eb = em.field.E(s2, pts);
tol = 1e-6 * max(em.vec.mag(Eb));
em.test.assertClose(max(abs(Eb(:,1))), 0, tol, 'bisector symmetry');

% 4. Dipole far-field scaling.
sd = em.src.merge(em.src.pointCharge(1e-9,[0 0 0.05]), ...
                  em.src.pointCharge(-1e-9,[0 0 -0.05]));
E1 = em.vec.mag(em.field.E(sd, [0 0 10]));
E2 = em.vec.mag(em.field.E(sd, [0 0 20]));
assert(abs(E1/E2 - 8)/8 <= 0.05, ...
    'Dipole far-field ratio must approach 8 within five percent.');

% 5. Merged charge element count is additive.
assert(size(s2.pos,1) == size(q1.pos,1) + size(q2.pos,1), ...
    'Merged element count must equal the sum of the inputs.');

% 6. Mixed charge/current sources must be rejected.
current = struct('type','current','pos',[0 0 0],'w',1, ...
    'Idl',[1 0 0],'label','test current');
caught = false;
try
    em.src.merge(q1, current);
catch
    caught = true;
end
assert(caught, 'Merging charge and current sources must raise an error.');

% 7. The field at a source location is singular and must raise an error.
caught = false;
try
    em.field.E(s, [0 0 0]);
catch
    caught = true;
end
assert(caught, 'Field evaluation at a point charge must raise an error.');

% 8. Shape and performance at N = 5000.
pts5k = randn(5000,3) + 5;
tic; E5k = em.field.E(s2, pts5k); elapsed = toc;
assert(isequal(size(E5k), [5000 3]), 'em.field.E must return Nx3.');
assert(elapsed < 2, 'em.field.E must evaluate 5000 points in under 2 s.');
fprintf('  5000-point field evaluation: %.4f s\n', elapsed);
disp('  test_lab03 checks complete');
