% test_lab08.m Dielectric fluxes, stored energy, capacitance and resistance.
run(fullfile(fileparts(mfilename('fullpath')),'lab08_studies.m'));
assert(abs(S.plate/S.plateRef-1)<.01,'parallel plate: d/w=.05');
assert(S.layerError<.01,'face-aligned constant layers: series formula');
assert(S.pH>1.8&&S.pH<2.2&&S.pA>.8&&S.pA<1.2,'graded two-layer harmonic/arithmetic orders');
assert(S.coaxErrors(end)<.02,'coax finest-grid error must be below 2 percent');
assert(all(diff(S.coaxErrors)<0)&&S.pCoax>.5&&S.pCoax<1.5,'staircased coax refinement');
assert(abs(S.coax100/S.coax(end)-1)<1e-10,'capacitance must be independent of voltage');
assert(abs(S.coaxOffset/S.coax(end)-1)<1e-10,'capacitance must be gauge independent');
assert(abs(S.coaxFilled/S.coax(end)-3)<1e-10,'uniform permittivity scaling');
assert(abs(S.bar/S.barRef-1)<1e-10,'uniform bar L/(sigma*w)');
assert(abs(S.RC/S.RCref-1)<1e-10,'uniform geometry RC=eps/sigma');
assert(S.tangentError<1e-10&&max(abs(S.EtBelow+1))<1e-10,'nonzero tangential E continuous');
assert(S.normalDspread<1e-10,'normal D conserved across dielectric interface');
assert(all(diff(S.Cw)>0)&&all(diff(S.Z0)<0),'microstrip capacitance/impedance trends');
assert(abs(S.widths(end)/S.widths(1)-10)<1e-12,'microstrip decade sweep');
assert(all(S.Cw>S.Cair)&&all(S.Z0>0),'dielectric loading and positive impedance');
assert(S.Cw(end)>S.ppWide,'wide strip fringe excess over parallel-plate estimate');
% Independent energy calculation for a linear field on a rectangular grid.
G=em.solve.grid2d([0 2],[0 3],17,12);W=em.derive.energyE(G,2*G.X-3*G.Y,5);
assert(abs(W/(.5*e0*5*13*6)-1)<1e-12,'energy integral and rectangular-grid weights');
% Voltage reversal leaves energy and capacitance invariant.
m=false(size(G.X));m(:,[1 end])=true;vf=zeros(size(m));vf(:,1)=2;
cp=em.derive.capacitance(G,2,m,vf,struct('method','direct'));
cn=em.derive.capacitance(G,2,m,-vf,struct('method','direct'));
assert(abs(cp/cn-1)<1e-12,'voltage reversal');
% All iterative choices must solve the variable-coefficient problem too.
g=em.solve.grid2d([0 1],[0 1],11,12);ep=1+3*(g.Y<.5);
mk=false(size(g.X));mk([1 end],:)=true;v=zeros(size(mk));v(end,:)=1;
vd=em.solve.poisson2d(g,0,ep,mk,v,struct('method','direct'));
for method={'jacobi','gs','sor'}
 vi=em.solve.poisson2d(g,0,ep,mk,v,struct('method',method{1},'tol',1e-11));
 assert(max(abs(vi(:)-vd(:)))<1e-8,'variable-coefficient iterative/direct agreement');
end
bad={@() em.derive.energyE(G,zeros(2),1), ...
 @() em.derive.energyE(G,G.X,-1), ...
 @() em.derive.capacitance(G,1,m,zeros(size(m))), ...
 @() em.derive.resistance(G,0,m,vf), ...
 @() em.solve.poisson2d(G,0,1,false(size(m)),vf), ...
 @() em.solve.poisson2d(G,0,1,m,vf,struct('faceMean','geometric'))};
for k=1:numel(bad)
 caught=false;try,bad{k}();catch,caught=true;end
 assert(caught,'invalid argument must raise an error');
end
fprintf('  Lab 8: plate error %.3e; constant-layer error %.3e\n',abs(S.plate/S.plateRef-1),S.layerError);
fprintf('  Graded-layer orders harmonic %.3f; arithmetic %.3f\n',S.pH,S.pA);
fprintf('  Coax finest error %.3f%%; order %.3f\n',100*S.coaxErrors(end),S.pCoax);
fprintf('  Voltage ratio %.10f; filled ratio %.10f\n',S.coax100/S.coax(end),S.coaxFilled/S.coax(end));
fprintf('  Bar R %.8f ohm m; tangential E jump %.3e V/m\n',S.bar,S.tangentError);
fprintf('  test_lab08 checks complete\n');
