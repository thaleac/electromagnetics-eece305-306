%% EECE 306 Lab 10 notebook. Ampere's Circuital Law, Curl, and the Magnetic Gauss Law
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore
%
% Run |publish('Lab10_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages. The coax cell below
% is the slowest computation of the semester so far, expect about a
% minute.

%% Setup
clear; close all;
addpath(fullfile(pwd,'tests'));

%% The coaxial cable as a bundle of filaments
% The API has no volume current constructor, and none is needed. A
% uniform volume current is a merge of many parallel line currents, each
% carrying its share of the total. The inner conductor carries +1 A, the
% shell returns it.
% Use the same calculations as the regression test.
S = lab10_studies();
fprintf('%d filaments inside, %d in the shell.\n', S.counts);

%% Circulation against enclosed current, region by region
rho = S.rho;
Circ = S.circulation;
IencRef = S.reference;
Ienc = IencRef(rho);
rFine = linspace(0,max(rho),401);
figure('Position',[100 100 620 330]);
plot(rho,Circ,'o',rFine,IencRef(rFine),'-'); grid on
xlabel('loop radius (m)'); ylabel('A')
legend('circulation of H','analytic enclosed current','Location','south')
title('Ampere''s law in all four regions of the coax')
fprintf('largest deviation = %.3f A\n',max(abs(Circ-Ienc)));
fprintf('largest outside circulation = %.3e A (consistent with zero)\n', ...
    max(abs(Circ(rho>S.edges(3)))));

%% One enclosed current, three path shapes
circle = S.paths(1);
sq = S.paths(2);
irr = S.paths(3);
fprintf('circle    %.5f A\nsquare    %.5f A\nirregular %.5f A\n',circle,sq,irr);

%% The magnetic Gauss law
% Compare the flux for a finite straight segment and a current loop.
% The hemisphere includes its closing disk. The reference scale is the
% field strength at the surface center times a disk or cube face area.
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

%% The point form
% Use the same points for curl H and div B. A 0.15 m step for the
% interior curl averages across the filament grid.
J = S.J;
cin = mean(S.curlIn,1);
cout = mean(S.curlOut,1);
divB = max(abs(reshape(S.divB,3,2)),[],1);
fprintf('Region                   Mean curl H (A/m^2)              Max |div B| (T/m)\n');
fprintf('Inner conductor       [%9.2e %9.2e %9.2e]     %12.3e\n',cin,divB(1));
fprintf('Dielectric            [%9.2e %9.2e %9.2e]     %12.3e\n',cout,divB(2));
fprintf('Reference Jz = %.4f A/m^2.\n',J);
%%
% The curl picked up the current in the conductor and dropped close to
% zero in the dielectric. The divergence stayed close to zero in both.

%% What field lines can and cannot show
wire1 = em.src.lineCurrent(1, @(t)[0.3+0*t 0*t t], [-40 40], 400);
wire2 = em.src.lineCurrent(-1, @(t)[-0.3+0*t 0*t t], [-40 40], 400);
pair = em.src.merge(wire1, wire2);
figure('Position',[100 100 580 340]);
em.viz.streamlines(@(r) em.field.H(wire1, r), [-1 1], [-1 1], ...
    [0.41 0; 0.51 0; 0.62 0; 0.73 0], struct('title', 'One wire, circles'));
figure('Position',[100 100 580 340]);
em.viz.streamlines(@(r) em.field.H(pair, r), [-1.2 1.2], [-1.2 1.2], ...
    [0.24 0; 0.185 0; 0.14 0; -0.24 0; -0.185 0; -0.14 0], ...
    struct('title', 'Antiparallel pair'));
%%
% In the first plot, the lines make circles around the wire. In the
% second plot, the lines still loop back, but the wires sit off center
% in the loops. The magnetic Gauss law says div B is zero, so there is
% no point where a magnetic field line starts or stops. In Lab 3, the
% electric displacement lines started on positive charges and ended on
% negative charges. We chose where to start tracing each line here,
% so closer lines on the page do not necessarily mean a stronger field.

%% Interpretation
% Moving the path outward in the inner conductor picked up more current.
% In the gap, the circle, square, and irregular path all gave 1 A because
% they enclosed the same current. The reading dropped as the path picked
% up more of the return current in the shield. Outside the shield, we
% measured less than a millionth of an amp. That was the clearest check
% for us because any mismatch between the two currents would show up there.

%% Problems encountered
% The readings near the conductor edges did not follow the smooth curve
% exactly. A small change in radius can pick up another row of filaments,
% so we checked the five percent limit farther from the edges. The other
% tricky case was the sphere that crossed the loop. The field changes
% quickly close to the wire, and the net flux got closer to zero when
% we used a finer grid.

%% Full test suite
testLog = evalc('testResults = runTests();');
testLines = strsplit(testLog,newline);
for k = 1:numel(testLines)
    if ~isempty(regexp(testLines{k},'^\s*(PASS|FAIL|\d+ passed)','once'))
        fprintf('%s\n',testLines{k});
    end
end
assert(~isempty(testResults) && all([testResults.passed]),'See the full runTests output for failures.');
