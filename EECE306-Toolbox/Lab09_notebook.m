%% EECE 306 Lab 9 notebook. Current Sources and the Biot-Savart Law
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore
%
% Current sources use the same element-list setup as the charge sources.
% The checks below compare the field with formulas, symmetry, and limits.

%% Setup
clear; close all;
addpath(pwd);
run(fullfile('tests','lab09_studies.m'));

%% How long em.field.H took
% em.field.H is the same loop as em.field.E with a different numerator.
minutes_for_H = 1.95;
fprintf('H review-to-first-passing-check interval: %.2f minutes\n',minutes_for_H);
% The recorded interval was about 1.9 minutes, from reviewing E through
% the first independent check of H, so it is not a timing for the whole lab.
% E already handled observation points in chunks and applied its constant
% at the end, so that part stayed the same. H reads Idl instead of q.*w,
% uses Idl cross (r-source), and divides by 4*pi instead of 4*pi*eps0.
% E itself did not need a rewrite; keeping the vector source data separate
% from the kernel made this much faster than starting over.

%% Circular loop against the closed form
a = 1; I = 1;
ring = @(t)[a*cos(t) a*sin(t) zeros(size(t))];
sr = em.src.lineCurrent(I, ring, [0 2*pi], 400);
z = linspace(0, 2, 10)';
Hnum = em.field.H(sr, [zeros(10,2) z]);
Href = I*a^2 ./ (2*(a^2 + z.^2).^1.5);
disp('    z (m)     relative error')
fprintf('%9.4f     %.3e\n',[z abs(Hnum(:,3)-Href)./Href]');
Hc = em.field.H(sr, [0 0 0]);
fprintf('center H_z = %.6f   I over 2a = %.6f\n', Hc(3), I/(2*a));
fprintf('transverse component at center = %.3e  (consistent with zero)\n', hypot(Hc(1), Hc(2)));

%% Finite segment against the infinite wire
% The same study as the Lab 4 finite line, with the field engine swapped.
% The curve has the same shape, and the handout asks why.
rho = 1;
hor = S9.Lrho; ratio = S9.wireRatio;
figure; semilogx(hor,ratio,'o-',hor,S9.chargeRatio,'x--'); grid on
xlabel('L / rho (total segment length)'); ylabel('finite / infinite field')
legend('current segment','Lab 4 line charge','Location','southeast')
title('Finite segment approaching the infinite wire')
fprintf('first sampled 99%%: L/rho %.3f; half-length h/rho %.3f\n', ...
    S9.wireThreshold,S9.wireThreshold/2);
fprintf('continuous threshold: L/rho %.3f; half-length h/rho %.3f\n', ...
    S9.wireCritical,S9.wireCritical/2);
% The magnetic and Lab 4 electric curves overlap because their transverse
% kernels have the same distance dependence. The starter called the full
% length h, while the handout uses 2h; here L is the full length, so the
% thresholds match only after that factor of two is accounted for.

%% Solenoid, built by reusing the ring
aS = S9.a; hS = S9.height; Nturns = S9.turns;
sol = S9.sol; zax = S9.zax; Hax = S9.Hax;
figure;
subplot(1,2,1); semilogx(S9.aspect,S9.solRatio,'o-'); grid on
xlabel('h / a'); ylabel('H center / (N I / h)')
title('Solenoid long-limit study')
subplot(1,2,2); plot(zax,Hax(:,3),'-'); grid on
xlabel('z (m)'); ylabel('H_z (A/m)')
title('Axial profile, 200 midpoint turns')
fprintf('solenoid first sampled 99%%: h/a %.1f\n',S9.solThreshold);
fprintf('center H_z %.5f; NI/h %.5f; ratio %.5f\n', ...
    S9.center,Nturns*I/hS,S9.center/(Nturns*I/hS));
fprintf('mouth / center %.5f; turn-refinement change %.3e\n', ...
    S9.mouthRatio,S9.turnRefinement);
% We used 200 evenly spaced midpoint turns instead of 25 endpoint turns.
% That keeps the turn density close to N/h and resolves the winding, so
% the 99% crossing is a long-solenoid limit rather than a spacing artifact.
% The mouth field is about half the center field, and it drops outside.

%% Reversal, linearity, and the merge policy
sneg = em.src.lineCurrent(-I, ring, [0 2*pi], 400);
fprintf('max |H(-I) + H(I)|  = %.3e  (reversal)\n', ...
    max(em.vec.mag(em.field.H(sneg, [0.3 0.2 0.5]) + em.field.H(sr, [0.3 0.2 0.5]))));
s2I = em.src.lineCurrent(2*I, ring, [0 2*pi], 400);
fprintf('max |H(2I) - 2H(I)| = %.3e  (linearity)\n', ...
    max(em.vec.mag(em.field.H(s2I, [0.3 0.2 0.5]) - 2*em.field.H(sr, [0.3 0.2 0.5]))));
sm = em.src.merge(sr, s2I);
fprintf('merged current source element count = %d  (sum of parts = %d)\n', ...
    size(sm.pos,1), size(sr.pos,1) + size(s2I.pos,1));
try
    em.src.merge(sr, em.src.pointCharge(1e-9, [0 0 0]));
    disp('ERROR, this line must not be reached')
catch err
    fprintf('merging current with charge correctly threw. %s\n', err.message);
end

%% Timing at N = 5000
pts = [5+sin((1:5000)'*.1), 5+cos((1:5000)'*.2), ...
    5+sin((1:5000)'*.3)];
tic; em.field.H(sr, pts); t = toc;
fprintf('em.field.H on 5000 points took %.3f s  (requirement, under 2 s)\n', t);

%% Interpretation
% The source list kept the geometry out of the field engine, and the
% quadrature code already knew how to turn a curve into weighted elements.
% We added an optional tangent output to quad.line without changing its
% existing positions or length weights, so lineCurrent reuses that work.
% The surface-current source uses K times area, while the line source keeps
% the tangent direction. For the next interface, we would keep the directed
% element data available instead of reducing everything to a magnitude.

%% Problems encountered
% The full-length versus half-length convention could have made the wire
% threshold look wrong by a factor of two. The first sampled crossing is
% only a bracket; it is not the exact continuous 99% aspect ratio.
% The original sparse endpoint-turn stack also changes the turn density,
% so we used midpoint turns and checked them against a finer winding.
% Near-zero transverse components need an absolute tolerance: a relative
% error there mostly measures division by cancellation noise.
% Reversing a current flips the arrows; B is mu0*mur times H, not H itself.

%% Full test suite
runTests
