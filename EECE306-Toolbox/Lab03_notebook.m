%% EECE 306 Lab 3 notebook. Point Charges, Superposition, and the Field Engine
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore
%
% Fill in every TODO, run |publish('Lab03_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages.

%% Setup
clear; close all;
rng(306);
e0 = em.const.eps0();

%% Closed form comparison at ten radii
Q = 2e-9;
s = em.src.pointCharge(Q, [0 0 0]);
radii = (1:10)';
Enum  = em.field.E(s, [radii zeros(10,2)]);
Eref  = Q ./ (4*pi*e0*radii.^2);
relerr = abs(Enum(:,1) - Eref) ./ Eref;
disp('   r (m)    relative error')
disp([radii relerr])

%% Symmetry on the bisector plane
% Two equal charges on the x axis. On the plane x = 0 the x component
% must vanish.
q1 = em.src.pointCharge(1e-9, [-0.1 0 0]);
q2 = em.src.pointCharge(1e-9, [ 0.1 0 0]);
s2 = em.src.merge(q1, q2);
pts = [zeros(200,1) randn(200,2)];
Ev  = em.field.E(s2, pts);
fprintf('max |Ex| on bisector      = %.3e V/m\n', max(abs(Ev(:,1))));
fprintf('max |E| on the same plane = %.3e V/m\n', max(em.vec.mag(Ev)));
% The expected x component is zero, so a relative tolerance cannot be used.
% We chose an absolute threshold of 1e-6 times the maximum field magnitude
% on the plane. This threshold is tiny compared with the field itself while
% still allowing for floating-point cancellation between the two charges.

%% Dipole far field
sd = em.src.merge(em.src.pointCharge(1e-9,[0 0 0.05]), em.src.pointCharge(-1e-9,[0 0 -0.05]));
E1 = em.vec.mag(em.field.E(sd, [0 0 10]));
E2 = em.vec.mag(em.field.E(sd, [0 0 20]));
fprintf('|E(r)| / |E(2r)| = %.4f   (expect near 8 for a dipole)\n', E1/E2);

%% The picture, raw and normalized
figure;
em.viz.quiver2(@(r) em.field.E(s2, r), [-0.5 0.5], [-0.5 0.5], 15, ...
    struct('normalize', false, 'title', 'Two equal charges, raw magnitudes'));
figure;
em.viz.quiver2(@(r) em.field.E(s2, r), [-0.5 0.5], [-0.5 0.5], 15, ...
    struct('normalize', true, 'title', 'Two equal charges, direction only'));
% The raw plot is difficult to read because the field strength changes by
% several orders of magnitude, so nearby arrows dominate the picture. The
% normalized plot clearly shows direction, but it gives up information about
% the relative field magnitude at each point.

%% Timing at N = 5000
pts5k = randn(5000,3) + 5;
tic; em.field.E(s2, pts5k); t = toc;
fprintf('em.field.E on 5000 points took %.3f s  (requirement, under 2 s)\n', t);

%% Interpretation
% The closed-form test checks the Coulomb constant, inverse-square magnitude,
% units, and radial components against a known analytical result. The
% symmetry test can reveal unequal treatment of mirrored elements or sign
% errors without needing a reference formula. The dipole scaling test checks
% whether superposition cancels the monopole term and produces the expected
% far-field 1/r^3 behavior. Each test therefore catches a different class of
% error that the other two could miss.

%% Problems encountered
% The supplied toolbox stored c0 without its required .m extension and did
% not include the cumulative Lab 2 test file. We restored those items, then
% corrected the Lab 3 implementation until the complete test suite passed.

%% Full test suite
runTests
