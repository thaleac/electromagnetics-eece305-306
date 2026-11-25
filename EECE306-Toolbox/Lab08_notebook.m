%% EECE 306 Lab 8 notebook. Dielectrics, Boundary Conditions, and Capacitance
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore
%
% This notebook collects the Lab 8 calculations, checks, and design curve.
% We used harmonic face permittivity for the solver and discrete face
% gradients for the stored energy, so both calculations use the same grid.

%% Setup
clear; close all;
addpath(pwd);
e0 = em.const.eps0();
run(fullfile('tests','lab08_studies.m'));

%% Rung 1. Parallel plate against the formula
% Plates span the full width, so the field between them is exactly one
% dimensional and C per unit depth must equal eps0 times w over d.
n = 82;
G = em.solve.grid2d([0 1], [0 0.1], n, n);
iy1 = round(0.25*(n-1)) + 1;  iy2 = round(0.75*(n-1)) + 1;
mask = false(n,n); mask(iy1,:) = true; mask(iy2,:) = true;
Vfix = zeros(n,n); Vfix(iy2,:) = 1;
C = em.derive.capacitance(G, 1, mask, Vfix, ...
    struct('method','direct'));
d = G.y(iy2)-G.y(iy1);
Cref = e0 / d;
fprintf('parallel plate C %.5e; formula %.5e; error %.2e\n', ...
    C, Cref, abs(C-Cref)/Cref);

%% Two dielectric layers in series
epsr = ones(n,n); epsr(G.Y < 0.05) = 4;
C2 = em.derive.capacitance(G, epsr, mask, Vfix, ...
    struct('method','direct'));
C2ref = e0 / ((d/2)/4 + (d/2)/1);
fprintf('two layers C %.5e; formula %.5e; error %.2e\n', ...
    C2, C2ref, abs(C2-C2ref)/C2ref);

% The plate gap is about 0.051 times the width. The side walls are insulating,
% so this is an ideal one-dimensional check with no edge fringing. It does
% not prove that a finite, open plate capacitor has negligible fringing.

%% Interface averaging and grid refinement
% The constant layers match the series result to round-off when the jump
% falls between nodes. There is no useful convergence slope left to fit.
% For the order check, each layer varies as (1+y), with a factor-of-four
% jump at y=0.5. The exact series integral is log(1.5)/4+log(2/1.5).
% This averaging comparison measures terminal face flux. The capacitance
% API itself uses stored energy with harmonic faces.
disp('N      harmonic relative error      arithmetic relative error')
fprintf('%4d       %12.4e              %12.4e\n',[S.N S.interfaceErrors]');
fprintf('graded-layer order: harmonic %.3f, arithmetic %.3f\n', S.pH, S.pA);
figure;
subplot(1,2,1); loglog(S.N-1,S.interfaceErrors,'o-'); grid on
xlabel('interval count'); ylabel('relative capacitance error')
legend('harmonic','arithmetic','Location','southwest')
title('Two graded dielectric layers')
subplot(1,2,2); loglog(S.coaxN-1,S.coaxErrors,'o-'); grid on
xlabel('interval count'); ylabel('relative capacitance error')
title('Staircased coaxial boundary')

%% Rung 2. Coaxial line against the log formula
n = 241;
Gc = em.solve.grid2d([-1 1], [-1 1], n, n);
R = sqrt(Gc.X.^2 + Gc.Y.^2);
a = 0.3; b = 0.9;
maskC = R <= a | R >= b;
VfixC = zeros(n,n); VfixC(R <= a) = 1;
Cc = em.derive.capacitance(Gc, 1, maskC, VfixC, ...
    struct('method','direct'));
CcRef = 2*pi*e0 / log(b/a);
fprintf('coax C %.5e; formula %.5e; error %.2e\n', ...
    Cc, CcRef, abs(Cc-CcRef)/CcRef);
disp('coax N      capacitance (F/m)      relative error')
fprintf('%4d       %12.5e           %12.4e\n',[S.coaxN S.coax S.coaxErrors]');
fprintf('coax observed order %.3f\n',S.pCoax);
% The coax capacitance comes up toward the formula as the grid gets finer.
% The circular boundaries are staircased, so the measured order is closer
% to one than two even though the interior stencil is second order.

%% The two blunder catchers
VfixC100 = zeros(n,n); VfixC100(R <= a) = 100;
Cc100 = em.derive.capacitance(Gc, 1, maskC, VfixC100, ...
    struct('method','direct'));
fprintf('C at 1 V: %.5e; C at 100 V: %.5e\n', Cc, Cc100);
% A missing voltage square would make these differ by 1e4.
Cfill = em.derive.capacitance(Gc, 3, maskC, VfixC, ...
    struct('method','direct'));
fprintf('epsr=3: C_filled/C = %.4f (expect 3)\n', Cfill/Cc);

%% Resistance of a uniform bar
n = 81;
Gr = em.solve.grid2d([0 1], [0 1], n, n);
maskR = false(n,n); maskR(:,1) = true; maskR(:,end) = true;
VfixR = zeros(n,n); VfixR(:,1) = 1;
sigma = 100;
Rbar = em.derive.resistance(Gr, sigma, maskR, VfixR, ...
    struct('method','direct'));
Rref = 1 / (sigma * 1);
fprintf('bar R %.5e ohm m; L/(sigma*w) %.5e ohm m\n', Rbar, Rref);

%% Rung 3. The payoff. A microstrip design curve
% The strip is h=0.2 m above ground, with epsr=4 below it and air above.
% The grid covers x=[-2,2], y=[0,2] m at 0.01 m spacing. The outer walls
% are grounded, so this is a boxed approximation to an open microstrip.
% We solve the same geometry with air to get C0, then use
% Z0 = 1/(c0*sqrt(C*C0)) for the lossless, nonmagnetic quasi-TEM model.
widths = S.widths; hsub = S.hsub; Cw = S.Cw;
figure;
subplot(1,2,1); plot(widths/hsub,Cw*1e12,'o-'); grid on
xlabel('w / h'); ylabel('C (pF/m)'); title('Microstrip capacitance')
subplot(1,2,2); plot(widths/hsub,S.Z0,'o-'); grid on
xlabel('w / h'); ylabel('Z0 (ohm)'); title('Quasi-TEM impedance')
disp('w/h       C (pF/m)       C0 air (pF/m)       Z0 (ohm)')
disp([widths(:)/hsub Cw(:)*1e12 S.Cair(:)*1e12 S.Z0(:)])
fprintf('widest strip: C %.3f pF/m; plate estimate %.3f pF/m\n', ...
    Cw(end)*1e12,S.ppWide*1e12);
fprintf('widest strip: mesh change %.2f%%; larger box change %.2f%%\n', ...
    100*S.meshChange,100*S.boxChange);
% Making the strip wider raises C and lowers Z0 over the decade we checked.
% At w/h=3, C is about 140.2 pF/m versus 106.3 pF/m from the plate estimate;
% the extra capacitance comes from fringing, with some influence from the box.
% Halving the grid spacing still changes the widest case by about 6.7%, so
% these are design estimates, not fully mesh-converged open-line values.
% Enlarging the box at fixed coarse spacing changes that case by about 0.22%.

%% Boundary conditions and what the numbers mean
fprintf('tangential E jump %.3e V/m; normal D relative spread %.3e\n', ...
    S.tangentError,S.normalDspread);
fprintf('uniform geometry RC %.6e s; eps/sigma %.6e s\n',S.RC,S.RCref);
% The tangential-field check uses a nonzero field on both sides of the
% interface, so a zero-field symmetry case cannot make it pass by accident.
% Normal D stays continuous without free interface charge. Conductivity
% replaces permittivity in the steady-current equation; for one uniform
% geometry, C per depth times R times depth is eps/sigma, in seconds.
% The energy sums harmonic face coefficients times squared voltage drops.
% Resistance instead sums current over faces leaving the high terminal.
% Geometry and interface placement still limit the accuracy, so a passing
% formula check does not make the boxed microstrip an exact open-line result.

%% Full regression suite, Lab 1 through Lab 8
runTests
