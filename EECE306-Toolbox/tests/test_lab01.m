% test_lab01.m  Lab 1 checks. Add to this file, do not delete existing checks.
% Run from the toolbox root with:  runTests

em.test.assertClose(em.const.c0(), 1/sqrt(em.const.mu0()*em.const.eps0()), 1e-12, 'c0 derived');
em.test.assertClose(em.const.eta0(), 376.730, 1e-3, 'eta0 about 377 ohm');
em.test.assertClose(em.vec.mag([3 4 0]), 5, 1e-12, 'mag 3-4-5');
em.test.assertClose(em.vec.unit([3 4 0]), [0.6 0.8 0], 1e-12, 'unit 3-4-5');
assert(all(isnan(em.vec.unit([0 0 0]))), 'unit of zero vector must be NaN');
em.test.assertClose(em.vec.angle([1 0 0],[0 1 0]), pi/2, 1e-12, 'angle 90 deg');
th = em.vec.angle([1 0 0],[-1 0 0]);
assert(isreal(th), 'angle must be real for antiparallel vectors');
em.test.assertClose(th, pi, 1e-12, 'angle 180 deg');
em.test.assertClose(em.vec.fromTo([1 1 1],[2 3 4]), [1 2 3], 1e-12, 'fromTo');

A = randn(100,3);
assert(isequal(size(em.vec.mag(A)),  [100 1]), 'mag must return Nx1');
assert(isequal(size(em.vec.unit(A)), [100 3]), 'unit must return Nx3');
assert(isequal(size(em.vec.angle(A, A)), [100 1]), 'angle must return Nx1');
assert(isequal(size(em.vec.fromTo(A, A)), [100 3]), 'fromTo must return Nx3');
disp('  test_lab01 checks complete');
