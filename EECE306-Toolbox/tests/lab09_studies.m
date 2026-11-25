% Reproducible source/field studies shared by Lab 9 notebook and tests.
S9=struct();a=1;I=1;ring=@(t)[a*cos(t) a*sin(t) zeros(size(t))];
S9.ring=em.src.lineCurrent(I,ring,[0 2*pi],400);
S9.z=linspace(0,2,10)';S9.axis=em.field.H(S9.ring,[zeros(10,2) S9.z]);
S9.ref=I*a^2./(2*(a^2+S9.z.^2).^1.5);
S9.loopError=max(abs(S9.axis(:,3)./S9.ref-1));
S9.transverse=max(max(abs(S9.axis(:,1:2))));
% Use total length L, matching Lab 4; the Lab 9 handout calls L/2 'h'.
S9.Lrho=logspace(0,2,15)';S9.wireRatio=zeros(15,1);S9.chargeRatio=zeros(15,1);
seg=@(t)[zeros(size(t)) zeros(size(t)) t];
for k=1:15
 L=S9.Lrho(k);s=em.src.lineCurrent(I,seg,[-L/2 L/2],800);
 H=em.field.H(s,[1 0 0]);S9.wireRatio(k)=H(2)/(I/(2*pi));
 q=em.src.lineCharge(1e-9,seg,[-L/2 L/2],800);
 E=em.field.E(q,[1 0 0]);S9.chargeRatio(k)=E(1)/(1e-9/(2*pi*em.const.eps0()));
end
S9.wireExact=S9.Lrho./sqrt(S9.Lrho.^2+4);
S9.wireFirst=find(S9.wireRatio>.99,1);
S9.wireThreshold=S9.Lrho(S9.wireFirst);
S9.wireCritical=2*.99/sqrt(1-.99^2);
% Finite winding, with midpoint turns and enough turns per radius to resolve it.
S9.aspect=[2 4 8 10 12 14 16 20 30 40]';S9.solRatio=zeros(10,1);
S9.turns=200;S9.a=.4;S9.height=8;phiN=80;
for k=1:10
 h=S9.a*S9.aspect(k);zk=(-.5+((1:S9.turns)-.5)/S9.turns)*h;
 srcs=cell(1,S9.turns);
 for j=1:S9.turns
  zc=zk(j);srcs{j}=em.src.lineCurrent(I,@(t)[S9.a*cos(t) S9.a*sin(t) zc*ones(size(t))],[0 2*pi],phiN);
 end
 s=em.src.merge(srcs{:});hc=em.field.H(s,[0 0 0]);
 S9.solRatio(k)=hc(3)/(S9.turns*I/h);
 if S9.aspect(k)==20,S9.sol=s;S9.zk=zk;end
end
S9.solFirst=find(S9.solRatio>.99,1);S9.solThreshold=S9.aspect(S9.solFirst);
S9.solExact=S9.aspect./sqrt(S9.aspect.^2+4);
S9.zax=linspace(-6,6,121)';S9.Hax=em.field.H(S9.sol,[zeros(121,2) S9.zax]);
S9.HaxisRef=zeros(size(S9.zax));
for j=1:S9.turns
 S9.HaxisRef=S9.HaxisRef+I*S9.a^2./(2*(S9.a^2+(S9.zax-S9.zk(j)).^2).^1.5);
end
center=em.field.H(S9.sol,[0 0 0]);mouth=em.field.H(S9.sol,[0 0 S9.height/2]);
S9.center=center(3);S9.mouthRatio=mouth(3)/center(3);
% Refining turn placement checks the sheet approximation independently.
zk=(-.5+((1:400)-.5)/400)*S9.height;
S9.fineCenter=sum(I*S9.a^2./(2*(S9.a^2+zk.^2).^1.5));
S9.turnRefinement=abs(S9.fineCenter/(2*S9.center)-1);
