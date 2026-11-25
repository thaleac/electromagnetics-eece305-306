%% EECE 306 Lab 10. Ampere's law, curl, and magnetic flux
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore

%% Ampere's law in a coaxial cable
% The inner conductor carries 1 A and the shield carries 1 A back in the
% opposite direction. We split each current among parallel filaments on a
% 0.10 m grid. This lets us use the same magnetic field calculation as in
% Lab 9 and the same path integral as in Lab 6.
% Run or publish this script from the EECE306-Toolbox root.
addpath(fullfile(pwd,'tests'));
S = lab10_studies();
r = linspace(0,max(S.rho),401);
figure('Position',[100 100 620 330]);
plot(r,S.reference(r),'-',S.rho,S.circulation,'o','LineWidth',1.2);
grid on
xlabel('Loop radius (m)'); ylabel('Current (A)');
legend('Analytic enclosed current','Circulation of H','Location','south');
title('Ampere''s law across the coax');
fprintf('Inner conductor: %d filaments. Shield: %d filaments.\n',S.counts);
fprintf('Largest sampled disagreement: %.3f A.\n',max(abs(S.circulation-S.reference(S.rho))));
fprintf('Largest circulation outside the shield: %.3e A.\n',max(abs(S.circulation(S.rho>1.3))));
fprintf('Circle: %.5f A   Square: %.5f A   Irregular path: %.5f A\n',S.paths);
%%
% Inside the core, the enclosed current rose with radius. It stayed near
% 1 A in the gap and fell as the loop picked up more return current in the
% shield. Outside the cable, the two currents canceled. That region is a
% good check because even a small imbalance would show up against the
% zero-current reference. The noncircular paths gave the same result as the
% circle because all three enclosed the same filaments.

%% Magnetic flux through closed surfaces
% We reused the Lab 5 flux integral. The hemisphere includes a flat disk,
% and each cube face points outward. The reference scale is the field
% strength at the surface center times the disk or face area.
labels = {'Sphere','Hemisphere + disk','Cube'};
sources = {'Straight segment','Current loop'};
fprintf('%-18s %-19s %12s %12s\n','Source','Surface','Flux (Wb)','Scale (Wb)');
for i = 1:2
    for j = 1:3
        fprintf('%-18s %-19s %12.3e %12.3e\n',sources{i},labels{j},S.flux(i,j),S.fluxScale(i,j));
    end
end
fprintf('%-18s %-19s %12.3e %12.3e\n','Current loop','Crossing sphere',S.cutFlux(end),S.cutScale);
fprintf('Crossing sphere, coarse and fine flux: %.3e and %.3e Wb.\n',S.cutFlux);
%%
% The crossing sphere cuts through the wire in two places. The field
% changes quickly nearby, so that row is more sensitive to the integration
% grid. Its net flux is still small relative to the reference scale.

%% The point form
% For the curl inside the core we used a 0.15 m step. A much smaller step
% would sample the gaps between filaments instead of averaging across them.
% The dielectric has no filaments, so we could use a smaller step there.
C = {S.curlIn,S.curlOut};
D = reshape(S.divB,3,2);
region = {'Inner conductor','Dielectric'};
fprintf('%-17s %10s %10s %12s %14s\n','Region','Curl Hx','Curl Hy','Mean curl Hz','Max |div B|');
for k = 1:2
    fprintf('%-17s %10.3e %10.3e %12.3e %14.3e\n',region{k},mean(C{k},1),max(abs(D(:,k))));
end
fprintf('Reference Jz = %.4f A/m^2. Curl is in A/m^2 and divergence in T/m.\n',S.J);
%%
% Curl picked up the current density in the inner conductor and dropped
% close to zero in the gap, while divergence stayed close to zero in both.

%% Field lines around one wire
wire1 = em.src.lineCurrent(1,@(t) [0.3+0*t,0*t,0.16*sinh(t)],[-6.2 6.2],240);
wire2 = em.src.lineCurrent(-1,@(t) [-0.3+0*t,0*t,0.16*sinh(t)],[-6.2 6.2],240);
figure('Position',[100 100 580 340]);
em.viz.streamlines(@(r) em.field.H(wire1,r),[-0.25 0.85],[-0.55 0.55], ...
    [0.41 0; 0.51 0; 0.62 0; 0.73 0],struct('title','Field around one wire'));

%% Field lines around two opposite currents
pair = em.src.merge(wire1,wire2);
figure('Position',[100 100 580 340]);
em.viz.streamlines(@(r) em.field.H(pair,r),[-1 1],[-0.65 0.65], ...
    [0.24 0; 0.185 0; 0.14 0; -0.24 0; -0.185 0; -0.14 0], ...
    struct('title','Field around two opposite currents'));
%%
% The paths around one wire are circles. With two wires the paths bend
% toward the gap and still loop back in the region shown. This fits the
% magnetic Gauss law: zero divergence means there are no sources or sinks
% where a magnetic field line can start or stop. In Lab 3, the electric
% displacement lines could start on a positive charge or end on a negative
% charge. The seed locations set the spacing in these plots, so the spacing
% alone cannot tell us the field strength.

%% Numerical issues and full test output
% The grid of filaments makes the enclosed current change in steps near a
% conductor edge. We kept those radii in the figure and used radii farther
% from the edges for the five percent comparison.
[testLog,testResults] = evalc('runTests');
testLines = strsplit(testLog,newline);
for k = 1:numel(testLines)
    if ~isempty(regexp(testLines{k},'^\s*(PASS|FAIL|\d+ passed)','once'))
        fprintf('%s\n',testLines{k});
    end
end
assert(~isempty(testResults) && all([testResults.passed]),'See the full runTests output for failures.');
