% Shared reproducible Lab 8 experiments. No tests are skipped or cached.
e0=em.const.eps0(); S=struct(); direct=struct('method','direct');
% Ideal full-width plates with insulating sides: exactly one-dimensional.
G=em.solve.grid2d([0 1],[0 .05],21,82); m=false(size(G.X));m([1 end],:)=true;
v=zeros(size(G.X));v(end,:)=1; e=ones(size(v));e(G.Y<.025)=4;
S.plate=em.derive.capacitance(G,1,m,v,direct);S.plateRef=e0/.05;
S.layer=em.derive.capacitance(G,e,m,v,direct);S.layerRef=e0/(.025/4+.025);
S.aspect=.05;S.layerError=abs(S.layer/S.layerRef-1);
% A smooth variation inside each layer leaves a measurable quadrature error.
% Constant layers with a face-aligned interface are exact to round-off.
S.N=[20 40 80 160]'; S.interfaceErrors=zeros(4,2);
for k=1:4
 n=S.N(k);G=em.solve.grid2d([0 1],[0 1],5,n);
 e=(1+G.Y).*(1+3*(G.Y<.5));m=false(n,5);m([1 end],:)=true;v=zeros(n,5);v(end,:)=1;
 for j=1:2
  kinds={'harmonic','arithmetic'};
  V=em.solve.poisson2d(G,0,e,m,v,struct('method','direct','faceMean',kinds{j}));
  if j==1,ef=2./(1./e(1,:)+1./e(2,:));else,ef=(e(1,:)+e(2,:))/2;end
  c=e0*sum(ef.*(V(2,:)-V(1,:)).*[.5 1 1 1 .5])*G.hx/G.hy;
  ref=e0/(log(1.5)/4+log(2/1.5));
  S.interfaceErrors(k,j)=abs(c/ref-1);
 end
end
% N-1 counts intervals, so the fit uses h=1/(N-1).
S.pH=em.test.convergence(@(N) S.interfaceErrors(S.N-1==N,1),S.N-1,struct('plot',false));
S.pA=em.test.convergence(@(N) S.interfaceErrors(S.N-1==N,2),S.N-1,struct('plot',false));
S.coaxN=[61 121 241]';S.coax=zeros(3,1);S.coaxRef=2*pi*e0/log(3);
for k=1:3
 n=S.coaxN(k);G=em.solve.grid2d([-1 1],[-1 1],n,n);r=hypot(G.X,G.Y);
 m=r<=.3|r>=.9;v=zeros(n);v(r<=.3)=1;
 S.coax(k)=em.derive.capacitance(G,1,m,v,direct);
end
S.coaxErrors=abs(S.coax/S.coaxRef-1);
S.pCoax=em.test.convergence(@(N) S.coaxErrors(S.coaxN-1==N),S.coaxN-1,struct('plot',false));
S.coax100=em.derive.capacitance(G,1,m,100*v,direct);
S.coaxFilled=em.derive.capacitance(G,3,m,v,direct);
S.coaxOffset=em.derive.capacitance(G,1,m,v+7,direct);
Gr=em.solve.grid2d([0 1],[0 1],41,41);mr=false(41);mr(:,[1 end])=true;
vr=zeros(41);vr(:,1)=1;S.bar=em.derive.resistance(Gr,100,mr,vr,direct);
S.barRef=.01;S.RC=S.bar*em.derive.capacitance(Gr,3,mr,vr,direct);
S.RCref=3*e0/100;
% Tangential E is nonzero; sampling a purely normal field would be vacuous.
Gt=em.solve.grid2d([0 1],[0 1],41,42);et=1+3*(Gt.Y<.5);
mt=false(size(Gt.X));mt([1 end],:)=true;mt(:,[1 end])=true;
vt=Gt.X;Vt=em.solve.poisson2d(Gt,0,et,mt,vt,direct);
S.EtBelow=-diff(Vt(21,:))/Gt.hx;S.EtAbove=-diff(Vt(22,:))/Gt.hx;
S.tangentError=max(abs(S.EtBelow-S.EtAbove));
% Normal D on consecutive faces is constant in the constant-layer problem.
Gn=em.solve.grid2d([0 1],[0 1],5,82);en=1+3*(Gn.Y<.5);
mn=false(82,5);mn([1 end],:)=true;vn=zeros(82,5);vn(end,:)=1;
Vn=em.solve.poisson2d(Gn,0,en,mn,vn,direct);
D=-e0*(2./(1./en(1:end-1,3)+1./en(2:end,3))).*diff(Vn(:,3))/Gn.hy;
S.normalDspread=(max(D)-min(D))/max(abs(D));
% Boxed microstrip approximation. The same box defines dielectric/air pairs.
S.widths=[.06 .1 .2 .3 .4 .6];S.hsub=.2;S.Cw=zeros(1,6);S.Cair=S.Cw;
Gm=em.solve.grid2d([-2 2],[0 2],401,201);epsm=1+3*(Gm.Y<.2);
[~,iy]=min(abs(Gm.y-.2));
for k=1:6
 mk=false(size(Gm.X));mk([1 end],:)=true;mk(:,[1 end])=true;
 strip=abs(Gm.x)<=S.widths(k)/2+1e-12;mk(iy,strip)=true;
 vf=zeros(size(mk));vf(iy,strip)=1;
 S.Cw(k)=em.derive.capacitance(Gm,epsm,mk,vf,direct);
 S.Cair(k)=em.derive.capacitance(Gm,1,mk,vf,direct);
end
S.Z0=1./(em.const.c0()*sqrt(S.Cw.*S.Cair));
S.ppWide=4*e0*S.widths(end)/S.hsub;
% Separate mesh and domain checks at the widest strip, at matched geometry.
S.sensitivity=zeros(1,2);
for j=1:2
 if j==1,gb=em.solve.grid2d([-2 2],[0 2],201,101);
 else,gb=em.solve.grid2d([-3 3],[0 3],301,151);end
 ep=1+3*(gb.Y<.2);[~,ii]=min(abs(gb.y-.2));
 mb=false(size(gb.X));mb([1 end],:)=true;mb(:,[1 end])=true;
 strips=abs(gb.x)<=.3+1e-12;mb(ii,strips)=true;
 vb=zeros(size(mb));vb(ii,strips)=1;
 S.sensitivity(j)=em.derive.capacitance(gb,ep,mb,vb,direct);
end
S.meshChange=abs(S.Cw(end)/S.sensitivity(1)-1);
S.boxChange=abs(S.sensitivity(2)/S.sensitivity(1)-1);
