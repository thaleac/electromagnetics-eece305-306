%% EECE 306 Lab 6 notebook. Electric Potential, the Gradient, and Path Independence
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore
%
% Fill in every TODO, run |publish('Lab06_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages.

%% Setup
clear; close all;
e0 = em.const.eps0();
Q  = 2e-9;
s  = em.src.pointCharge(Q, [0 0 0]);
FE = @(r) em.field.E(s, r);
FV = @(r) em.field.V(s, r);

%% Potential against the closed form
radii = (1:10)';
Vnum = em.field.V(s, [radii zeros(10,2)]);
Vref = Q ./ (4*pi*e0*radii);
disp('   r (m)    relative error')
disp([radii abs(Vnum - Vref)./Vref])

%% E from the gradient, and the error map
% E computed directly and E recovered as minus the gradient of V must
% agree away from the source.
[X, Y] = meshgrid(linspace(0.3, 2, 40), linspace(-1, 1, 40));
P = [X(:) Y(:) zeros(numel(X),1)];
sMap = em.src.merge( ...
    em.src.pointCharge( Q, [-0.5 0 0]), ...
    em.src.pointCharge(-Q, [ 0.5 0 0]), ...
    em.src.pointCharge(Q/2, [0 0.5 0]));
FVMap = @(r) em.field.V(sMap, r);
Ed = em.field.E(sMap, P);
Eg = -em.op.grad(FVMap, P);
absoluteError = em.vec.mag(Ed - Eg);
fieldMagnitude = em.vec.mag(Ed);
relativeError = absoluteError ./ max(fieldMagnitude, realmin);
errmap = reshape(max(relativeError, realmin), size(X));
figure; contourf(X, Y, log10(errmap), 20); colorbar
xlabel('x (m)'); ylabel('y (m)')
title('log10 relative error of E recovered from  -grad V')
% The largest absolute discrepancy is 7.776e-4 V/m near (0.518,-0.026,0) m,
% close to the negative charge where rapid variation increases truncation error.
% The retained relative map shows the same location; the absolute value is
% reported separately because the handout asks for error magnitude.

%% Choosing the step size, the V shaped curve
% Truncation error falls as h squared, round off grows as 1 over h. The
% total has a minimum. Sweep h and find it for this machine.
p0 = [1 0.5 0];
Eexact = em.field.E(s, p0);
hs = logspace(-12, -1, 23)';
errh = zeros(size(hs));
for k = 1:numel(hs)
    Eh = -em.op.grad(FV, p0, hs(k));
    errh(k) = em.vec.mag(Eh - Eexact) / em.vec.mag(Eexact);
end
figure; loglog(hs, errh, 'o-'); grid on
xlabel('step h (m)'); ylabel('relative error of  -grad V')
title('The two error sources compete, the total has a minimum')
[emin, imin] = min(errh);
fprintf('best h on this machine about %.1e with error %.1e\n', hs(imin), emin);
fprintf('the API default h = 1e-5 sits on the flat bottom of this curve\n');

%% The field is conservative, with a control
qSquare = @(t) 2*t/pi;
square = @(t)[ ...
    -1+2*min(qSquare(t),1)-2*min(max(qSquare(t)-2,0),1), ...
    -1+2*min(max(qSquare(t)-1,0),1)-2*min(max(qSquare(t)-3,0),1), ...
    zeros(size(t))];
loops = { @(t)[2*cos(t) 2*sin(t) zeros(size(t))], square, ...
          @(t)[1.5+0.5*cos(t) 0.5*sin(3*t) 0.3*sin(t)] };
pathNames = {'circle', 'square', 'nonplanar loop'};
lineC = @(t)[zeros(size(t)) zeros(size(t)) t];
sLine = em.src.lineCharge(1e-9, lineC, [-0.3 0.3], 200);
ball = @(u,v,w)[w*sin(u)*cos(v) w*sin(u)*sin(v) w*cos(u)];
sVol = em.src.volCharge(1e-9, ball, [0 pi], [0 2*pi], [0 0.25], 8, 16, 8);
sources = {sMap, sLine, sVol};
sourceNames = {'point charges', 'finite line', 'volume charge'};
for j = 1:3
    Fsource = @(r) em.field.E(sources{j}, r);
    for k = 1:3
        C = em.field.circulation(Fsource, loops{k}, [0 2*pi], 400);
        fprintf('circulation of E: %s / %s = %.3e V\n', ...
            sourceNames{j}, pathNames{k}, C);
    end
end
Fnc = @(r) [-r(:,2) r(:,1) zeros(size(r,1),1)];
Cnc = em.field.circulation(Fnc, loops{1}, [0 2*pi], 400);
fprintf('circulation of the control field = %.4f  (expect 8 pi, clearly nonzero)\n', Cnc);
% The control matters. A test that only ever confirms zero cannot tell
% correct code from code that always returns zero.

%% Potential difference two ways
pa = [1 0 0]; pb = [0.3 0.4 0];
dV = em.field.V(s, pb) - em.field.V(s, pa);
seg = @(t) pa + t.*(pb - pa);
work = -em.field.circulation(FE, seg, [0 1], 400);
fprintf('V(b) - V(a) from em.field.V        = %.6e V\n', dV);
fprintf('minus line integral of E, a to b   = %.6e V\n', work);

%% Seeing the potential
figure;
sDip = em.src.merge( ...
    em.src.pointCharge( Q, [-0.5 0 0]), ...
    em.src.pointCharge(-Q, [ 0.5 0 0]));
em.viz.contourV(@(r) em.field.V(sDip, r), [-2 2], [-2 2], 60, ...
    struct('title', 'Equipotentials of a dipole'));
hold on
em.viz.quiver2(@(r) em.field.E(sDip, r), [-2 2], [-2 2], 16, ...
    struct('normalize', true, 'title', 'Equipotentials with field lines overlaid'));
title('Equipotentials with field lines overlaid')

%% Interpretation
% The sampled minimum is h=3.162e-6 m with relative error 1.025e-11;
% h=1e-5 m gives 3.586e-11, so it is a useful nearby compromise.
% Every later gradient should be checked at its own spatial scale because
% small h eventually amplifies subtraction round-off instead of improving accuracy.
% All nine closed integrals are below 3.5e-11 V, while the control gives 8*pi.
% The dipole arrows point toward decreasing potential; a sign error would
% reverse their direction toward the positive charge while keeping right angles.

%% Problems encountered
% The starter states the default is on a flat bottom, but the sampled minimum
% is at 3.162e-6 m and the default error is 3.50 times larger.
% The relative map is preserved; realmin guards zero magnitudes and log10(0),
% while tests measure the unmodified absolute discrepancy.
% A finite element list cannot identify infinite-source intent, so V warns
% when optional logical s.extentGrowsWithN is true (set after merging).
% Fixed-geometry refinement changes V by 4.965e-6 V without a warning;
% doubling line extent changes V from 25.9495 to 37.6527 V.
% The 400-panel open integral differs by 2.706e-5 V; doubling panels
% quarters this error, and an independent Gauss integral agrees within 1e-8 V.
% Square corners align with midpoint panels; exact-tangent edge sums
% agree within 4.2e-11 V, avoiding an assumption of smoothness at corners.

%% Full test suite
runTests
