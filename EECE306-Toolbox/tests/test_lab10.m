% Lab 10 regression checks for Ampere's law, flux, curl and field lines.
S10 = lab10_studies();
rho10 = S10.rho;
rEdges10 = S10.edges;
away10 = all(abs(rho10-rEdges10) >= S10.spacing-1e-12,2) & rho10>0;
region10 = 1+sum(rho10>rEdges10,2);
assert(isequal(unique(region10(away10)),(1:4)'), 'Expected samples in all four coax regions.');
assert(all(abs(S10.circulation(away10)-S10.reference(rho10(away10))) <= 0.05*S10.I), ...
    'Coax circulation disagrees with enclosed current.');
% Equal and opposite filaments nearly cancel outside. Allow 0.01% of I
% for the finite-length end field and numerical quadrature.
assert(all(abs(S10.circulation(region10==4)) < 1e-4*S10.I), 'Outside coax circulation.');
assert(max(S10.paths)-min(S10.paths) <= 0.01*S10.I, 'Path shape changes circulation.');
assert(all(abs(S10.paths-S10.I) <= 0.01*S10.I), 'Paths should enclose the full inner current.');

assert(all(abs(S10.flux(:)) < 1e-6*S10.fluxScale(:)), 'Closed-surface magnetic flux.');
assert(abs(S10.cutFlux(end)) < 0.005*S10.cutScale, 'Flux through surface crossing the loop.');

assert(all(em.vec.mag(S10.curlIn-[0 0 S10.J]) < 0.05*S10.J), 'Interior curl of H differs from J.');
assert(all(em.vec.mag(S10.curlOut) < 0.02*S10.J), 'Dielectric curl of H should vanish.');
assert(all(abs(S10.divB) < 1e-5*em.const.mu0()*S10.J), 'Divergence of coax B should vanish.');
assert(all(em.vec.mag(S10.curlE) < 1e-6), 'Electrostatic curl should vanish.');

% Different component functions catch a swapped sign or cyclic index.
positions10 = [1 2 3; -2 0 1];
F10 = @(r) [2*r(:,2)+3*r(:,3),5*r(:,1)+7*r(:,3),11*r(:,1)+13*r(:,2)];
em.test.assertClose(em.op.curl(F10,positions10),repmat([6 -8 3],2,1),1e-9,'Cartesian curl');
assert(isequal(size(em.op.curl(F10,zeros(0,3))),[0 3]));

fig10 = figure('Visible','off');
cleanup10 = onCleanup(@() close(fig10));
line10 = em.viz.streamlines(@(r) [-r(:,2),r(:,1),0*r(:,3)],[-1 1],[-1 1],[0.5 0]);
x10 = get(line10,'XData'); y10 = get(line10,'YData');
assert(numel(x10)>100 && all(abs(hypot(x10,y10)-0.5)<1e-4));
assert(hypot(x10(end)-x10(1),y10(end)-y10(1))<1e-9,'Circular line should close.');
line10 = em.viz.streamlines(@(r) [1+0*r(:,1),0*r(:,2),0*r(:,3)],[-1 1],[-1 1],[0 0]);
x10 = get(line10,'XData');
assert(min(x10)<-0.99 && max(x10)>0.99,'Open line should reach both plot boundaries.');
clear cleanup10
fprintf('  Lab 10: outside circulation %.2e A, crossing flux %.2e of reference scale\n', ...
    max(abs(S10.circulation(region10==4))),abs(S10.cutFlux(end))/S10.cutScale);
