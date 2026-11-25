% test_lab02.m  Lab 2 cumulative coordinate-system checks.
rng(306);
r = randn(1000,3);
A = randn(1000,3);

em.test.assertClose(max(em.vec.mag(em.coord.sph2c(em.coord.c2sph(r))-r)), ...
    0, 1e-10, 'spherical point roundtrip');
em.test.assertClose(max(em.vec.mag(em.coord.cyl2c(em.coord.c2cyl(r))-r)), ...
    0, 1e-10, 'cylindrical point roundtrip');
em.test.assertClose(max(em.vec.mag(em.coord.vecSph2C( ...
    em.coord.vecC2Sph(A,r),r)-A)), 0, 1e-10, 'spherical vector roundtrip');
em.test.assertClose(max(em.vec.mag(em.coord.vecCyl2C( ...
    em.coord.vecC2Cyl(A,r),r)-A)), 0, 1e-10, 'cylindrical vector roundtrip');

em.test.assertClose(em.coord.c2sph([0 0 1]), [1 0 0], 1e-12, 'north pole');
em.test.assertClose(em.coord.c2sph([0 0 -1]), [1 pi 0], 1e-12, 'south pole');
q = [1 1 0; -1 1 0; -1 -1 0; 1 -1 0];
p = em.coord.c2cyl(q);
assert(all(p(:,2) >= 0 & p(:,2) < 2*pi), 'phi must lie in [0,2*pi)');

As = em.coord.vecC2Sph(A,r);
Ac = em.coord.vecC2Cyl(A,r);
em.test.assertClose(max(abs(em.vec.mag(As)-em.vec.mag(A))), ...
    0, 1e-10, 'spherical magnitude invariant');
em.test.assertClose(max(abs(em.vec.mag(Ac)-em.vec.mag(A))), ...
    0, 1e-10, 'cylindrical magnitude invariant');
em.test.assertClose(em.coord.vecC2Sph([2 3 0],[0 5 0]), ...
    [3 0 -2], 1e-12, 'two-position case one');
em.test.assertClose(em.coord.vecC2Sph([2 3 0],[4 0 0]), ...
    [2 0 3], 1e-12, 'two-position case two');

% Off-axis reference case catches an incorrect spherical azimuth.
em.test.assertClose(em.coord.vecC2Sph([7 -2 5],[3 4 12]), ...
    [73/13 31/65 -34/5], 1e-12, 'off-axis spherical components');
disp('  test_lab02 checks complete');
