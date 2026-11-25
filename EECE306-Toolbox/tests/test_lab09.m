% test_lab09.m Current-source contracts, Biot-Savart signs, and closed forms.
run(fullfile(fileparts(mfilename('fullpath')),'lab09_studies.m'));
assert(S9.loopError<1e-6,'loop axis closed form at ten z values');
em.test.assertClose(S9.axis(1,3),.5,1e-6,'loop center I/(2a)');
% 1e-12 A/m is 2e-12 of the largest on-axis field (0.5 A/m),
% allowing cancellation/derivative round-off without using relative error at zero.
assert(S9.transverse<1e-12,'loop transverse symmetry in A/m');
assert(S9.wireRatio(S9.wireFirst)>.99&&S9.wireRatio(S9.wireFirst-1)<.99,'reported wire threshold bracket');
assert(max(abs(S9.wireRatio-S9.wireExact))<1e-6,'finite segment closed form');
assert(max(abs(S9.wireRatio-S9.chargeRatio))<1e-11,'Lab 4 and 9 dimensionless wire curves');
assert(S9.solRatio(S9.solFirst)>.99&&S9.solRatio(S9.solFirst-1)<.99,'reported solenoid threshold bracket');
assert(max(abs(S9.solRatio-S9.solExact))<1e-5,'dense turns approach finite-sheet formula');
assert(max(abs(S9.Hax(:,3)-S9.HaxisRef))<1e-8,'solenoid profile independent sum of closed-form loop fields');
assert(abs(S9.mouthRatio-.5)<.005&&S9.turnRefinement<1e-5,'mouth field and turn refinement');
ring=@(t)[cos(t) sin(t) zeros(size(t))];pts=[.3 .2 .5;-.5 .1 .7;.2 .4 -1];
s=S9.ring;h=em.field.H(s,pts);
sn=em.src.lineCurrent(-1,ring,[0 2*pi],400);s2=em.src.lineCurrent(2,ring,[0 2*pi],400);
em.test.assertClose(em.field.H(sn,pts),-h,1e-12,'componentwise current reversal');
em.test.assertClose(em.field.H(s2,pts),2*h,1e-12,'componentwise linearity');
sr=em.src.lineCurrent(1,@(t)ring(2*pi-t),[0 2*pi],400);
hr=em.field.H(sr,pts);
assert(max(abs(hr(:)+h(:)))<1e-10,'curve direction reversal with absolute tolerance near zero');
sm=em.src.merge(s,s2);assert(size(sm.pos,1)==800&&isequal(size(sm.Idl),[800 3]),'current merge element shapes');
em.test.assertClose(em.field.H(sm,pts),3*h,1e-12,'current merge superposition');
% Straight directed line: original two outputs and new tangent agree.
for rule={'midpoint','trapz','simpson','gauss'}
 c=@(t)[t 2*t 0*t];[p,w,tan]=em.quad.line(c,[0 2],21,rule{1});
 sl=em.src.lineCurrent(3,c,[0 2],21,rule{1});
 assert(max(abs(tan(:,1)-1))<1e-9&&max(abs(tan(:,2)-2))<1e-9,'raw tangent');
 em.test.assertClose(sl.pos,p,1e-12,'shared line positions');
 em.test.assertClose(sl.w,w,1e-12,'unchanged line weights');
 em.test.assertClose(sum(sl.Idl,1),[6 12 0],1e-9,'integrated I dr with all quadratures');
end
flat=em.src.lineCurrent(1,@(t)[0*t 0*t 0*t],[0 1],3);
assert(all(flat.Idl(:)==0)&&all(flat.w==0),'stationary curve must not divide by zero');
% Surface current: exact vector area, nonuniform and scalar-only K handles.
plane=@(u,v)[u v 0*u];sheet=em.src.surfCurrent(@(r)[2 0 0],plane,[-1 1],[-1 1],30,30);
em.test.assertClose(sum(sheet.Idl,1),[8 0 0],1e-9,'surface K times area');
nonuniform=em.src.surfCurrent(@(r)[r(1) 2*r(2) 0],plane,[-1 1],[-1 1],5,7);
em.test.assertClose(nonuniform.Idl(:,1),nonuniform.pos(:,1).*nonuniform.w,1e-12,'scalar-only K at each element');
vectorK=em.src.surfCurrent(@(r)[r(:,1) 2*r(:,2) zeros(size(r,1),1)],plane,[-1 1],[-1 1],5,7);
assert(max(abs(vectorK.Idl(:)-nonuniform.Idl(:)))<1e-12,'vectorized and scalar-only surface K agree');
flipped=em.src.surfCurrent(@(r)[2 0 0],@(u,v)[-u v 0*u],[-1 1],[-1 1],30,30);
hf=em.field.H(flipped,[.2 .3 1]);hp=em.field.H(sheet,[.2 .3 1]);
assert(max(abs(hf-hp))<1e-11,'surface coordinate reversal leaves physical K unchanged');
% Independent solid-angle field of a finite square current sheet at z=1.
Hsheet=em.field.H(sheet,[0 0 1]);theta=atan(1/sqrt(3));
em.test.assertClose(Hsheet,[0 -2*theta/pi 0],.001,'finite square sheet field direction and magnitude');
% Off-axis quadrature refinement rather than an axial symmetry-only test.
Ns=[20 40 80];er=zeros(1,3);ref=em.field.H(em.src.lineCurrent(1,ring,[0 2*pi],1280),[.9 0 .15]);
for k=1:3
 hh=em.field.H(em.src.lineCurrent(1,ring,[0 2*pi],Ns(k)),[.9 0 .15]);er(k)=norm(hh-ref);
end
assert(all(diff(er)<0)&&er(end)<er(1)/100,'off-axis loop refinement');
% Check the cross-product kernel with unrelated source vectors and weights.
st=struct('type','current','pos',[0 0 0;1 0 0],'w',[17;3],'Idl',[1 2 3;-2 1 4]);
p=[2 3 4;-1 2 3];reference=zeros(2,3);
for k=1:2
 d=p-st.pos(k,:);reference=reference+cross(repmat(st.Idl(k,:),2,1),d,2)./sum(d.^2,2).^1.5/(4*pi);
end
em.test.assertClose(em.field.H(st,p),reference,1e-12,'independent vector sum, no extra w');
assert(isequal(size(em.field.H(s,zeros(0,3))),[0 3]),'empty field dimensions');
em.test.assertClose(em.field.B(s,pts),em.const.mu0()*h,1e-12,'B default mur');
em.test.assertClose(em.field.B(s,pts,4),4*em.const.mu0()*h,1e-12,'B mur scaling');
bad={@()em.src.merge(s,em.src.pointCharge(1e-9,[0 0 0])), ...
 @()em.field.H(s,[1 2]),@()em.field.H(em.src.pointCharge(1,[0 0 0]),pts), ...
 @()em.field.H(s,s.pos(2,:)),@()em.field.B(s,pts,-1), ...
 @()em.src.lineCurrent([1 2],ring,[0 1],3), ...
 @()em.src.surfCurrent(@(r)1,plane,[0 1],[0 1],3,3)};
for k=1:numel(bad)
 caught=false;try,bad{k}();catch err
 caught=true;if k==4,assert(~isempty(strfind(err.message,'element 2')),'singularity must name element');end
 end
 assert(caught,'invalid arguments or mixed types must raise error');
end
pp=[5+sin((1:10000)'*.1),5+cos((1:10000)'*.2),5+sin((1:10000)'*.3)];
tic;hh=em.field.H(s,pp(1:5000,:));t5000=toc;
assert(isequal(size(hh),[5000 3])&&all(isfinite(hh(:)))&&t5000<2,'5000 positions in under 2 seconds');
tic;hh=em.field.H(s,pp);t10000=toc;
assert(isequal(size(hh),[10000 3])&&t10000<2,'10000 vectorized positions in under 2 seconds');
fprintf('  Lab 9: loop max relative error %.3e; transverse %.3e A/m\n',S9.loopError,S9.transverse);
fprintf('  Wire first 99%%: L/rho %.3f (half-length h/rho %.3f)\n',S9.wireThreshold,S9.wireThreshold/2);
fprintf('  Solenoid first 99%%: h/a %.1f; mouth/center %.6f\n',S9.solThreshold,S9.mouthRatio);
fprintf('  H timing: 5000 points %.3f s; 10000 points %.3f s\n',t5000,t10000);
fprintf('  test_lab09 checks complete\n');
