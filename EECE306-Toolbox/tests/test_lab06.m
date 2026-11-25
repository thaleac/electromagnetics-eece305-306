% test_lab06.m  Potential, derivatives, directed paths, and reference checks.
Q=2e-9; e0=em.const.eps0(); s=em.src.pointCharge(Q,[0 0 0]);
radii=(1:10)'; pts=[radii zeros(10,2)];
em.test.assertClose(em.field.V(s,pts),Q./(4*pi*e0*radii),2e-15,'ten radii');
assert(isequal(size(em.field.V(s,zeros(0,3))),[0 1]));
sm=em.src.merge(em.src.pointCharge(Q,[-0.5 0 0]), ...
    em.src.pointCharge(-Q,[0.5 0 0]),em.src.pointCharge(Q/2,[0 0.5 0]));
f=@(r) em.field.V(sm,r); p=[1 0.3 0.2; -1 0.5 0.3; 0.2 -1 0.4];
Ed=em.field.E(sm,p); Eg=-em.op.grad(f,p);
assert(isequal(size(Eg),[3 3]));
assert(max(em.vec.mag(Ed-Eg)./em.vec.mag(Ed))<1e-8);
vp=f(p);
for k=1:3, em.test.assertClose(em.field.V(sm,p(k,:)),vp(k),1e-14,'vectorized V'); end
% Polynomial derivatives independently expose component, sign, and 2h bugs.
fpoly=@(r) r(:,1).^3+2*r(:,2).^3+3*r(:,3).^3+7*r(:,1).*r(:,2);
gexact=[3*p(:,1).^2+7*p(:,2),6*p(:,2).^2+7*p(:,1),9*p(:,3).^2];
a=em.op.grad(fpoly,p,0.02)-gexact; b=em.op.grad(fpoly,p,0.01)-gexact;
em.test.assertClose(a,repmat([1 2 3]*0.02^2,3,1),1e-7,'central truncation');
em.test.assertClose(a./b,4,1e-6,'second order');
assert(max(abs(em.op.grad(fpoly,p)-gexact),[],'all')<2e-8);
assert(isequal(size(em.op.grad(fpoly,zeros(0,3))),[0 3]));
q=@(t) 2*t/pi;
square=@(t) [-1+2*min(q(t),1)-2*min(max(q(t)-2,0),1), ...
    -1+2*min(max(q(t)-1,0),1)-2*min(max(q(t)-3,0),1),zeros(size(t))];
loops={@(t) [2*cos(t) 2*sin(t) zeros(size(t))],square, ...
    @(t) [1.5+0.5*cos(t) 0.5*sin(3*t) 0.3*sin(t)]};
vertices=[-1 -1 0;1 -1 0;1 1 0;-1 1 0;-1 -1 0];
assert(max(abs(square((0:4)'*pi/2)-vertices),[],'all')<1e-14);
[~,ww]=em.quad.line(square,[0 2*pi],400);
assert(abs(sum(ww)-8)<1e-8);
for k=1:3, assert(norm(loops{k}(0)-loops{k}(2*pi))<1e-14); end
cloud=loops{3}(linspace(0,2*pi,100)'); assert(rank(cloud-mean(cloud))==3);
lc=@(t) [zeros(size(t)) zeros(size(t)) t];
sl=em.src.lineCharge(1e-9,lc,[-0.3 0.3],200);
ball=@(u,v,w) [w*sin(u)*cos(v) w*sin(u)*sin(v) w*cos(u)];
sv=em.src.volCharge(1e-9,ball,[0 pi],[0 2*pi],[0 0.25],8,16,8);
sources={sm,sl,sv}; ratios=zeros(3,3); squareDifference=0;
for j=1:3
    FF=@(r) em.field.E(sources{j},r);
    for k=1:3
        residual=zeros(1,3);
        for m=1:3
            NN=200*2^(m-1); residual(m)=em.field.circulation(FF,loops{k},[0 2*pi],NN);
        end
        [rp,wp]=em.quad.line(loops{k},[0 2*pi],800);
        fieldLengthScale=max(em.vec.mag(FF(rp)))*sum(wp);
        tol=1e-6*fieldLengthScale;
        assert(all(abs(residual)<tol),'closed path absolute field-times-length tolerance');
        assert(abs(residual(3)-residual(2))<tol/4,'closed path refinement');
        ratios(j,k)=abs(residual(3))/fieldLengthScale;
        reverse=@(t) loops{k}(2*pi-t);
        assert(abs(em.field.circulation(FF,reverse,[0 2*pi],800)+residual(3))<tol/4);
    end
    % Independent edge integral: exact displacements, no numerical tangent.
    edgeSum=0; [tt,wt]=em.quad.nodes(0,1,200,'midpoint');
    for edge=1:4
        delta=vertices(edge+1,:)-vertices(edge,:);
        edgeSum=edgeSum+sum(sum(FF(vertices(edge,:)+tt.*delta).*delta,2).*wt);
    end
    sqIntegral=em.field.circulation(FF,square,[0 2*pi],800);
    squareDifference=max(squareDifference,abs(sqIntegral-edgeSum));
    assert(abs(sqIntegral-edgeSum)<1e-8);
end
Fnc=@(r) [-r(:,2) r(:,1) zeros(size(r,1),1)];
for rule={'midpoint','trapz','simpson','gauss'}
    cn=em.field.circulation(Fnc,loops{1},[0 2*pi],80,rule{1});
    assert(abs(cn-8*pi)<2e-8,'nonzero 8pi control');
end
assert(abs(em.field.circulation(Fnc,@(t) loops{1}(2*pi-t),[0 2*pi],400)+8*pi)<2e-8);
assert(abs(em.field.circulation(Fnc,square,[0 2*pi],400)-8)<1e-8);
pa=[1 0 0]; pb=[0.3 0.4 0]; seg=@(t) pa+t.*(pb-pa);
dV=em.field.V(s,pb)-em.field.V(s,pa); FE=@(r) em.field.E(s,r);
Ns=[100 200 400 800]; er=zeros(size(Ns));
for k=1:4, er(k)=abs(-em.field.circulation(FE,seg,[0 1],Ns(k))-dV); end
assert(all(er(1:3)./er(2:4)>3.9) && all(er(1:3)./er(2:4)<4.1));
assert(er(end)/abs(dV)<1e-6);
assert(abs(-em.field.circulation(FE,seg,[0 1],32,'gauss')-dV)<1e-8);
% Fixed geometry refinement must not warn, regardless of count or amplitude.
lastwarn(''); vfinite=zeros(1,2); vgrow=zeros(1,2);
for k=1:2
    N=400*2^(k-1);
    fixed=em.src.lineCharge(1e-9,lc,[-1 1],N);
    vfinite(k)=em.field.V(fixed,[1 0 0]);
    growing=em.src.lineCharge(1e-9,lc,[-N/200 N/200],N);
    vgrow(k)=em.field.V(growing,[1 0 0]);
end
[~,wid]=lastwarn; assert(isempty(wid));
assert(abs(diff(vfinite))<1e-4 && diff(vgrow)>10);
growing.extentGrowsWithN=true;
warnState=warning('query','em:field:V:InfiniteReference');
warning('on','em:field:V:InfiniteReference'); lastwarn('');
warningText=evalc('em.field.V(growing,[1 0 0]);');
[~,wid]=lastwarn; warning(warnState.state,'em:field:V:InfiniteReference');
assert(strcmp(wid,'em:field:V:InfiniteReference'));
% Error contracts and source index propagation.
actions={@() em.field.V(sm,[0.5 0 0]), ...
    @() em.field.circulation(@(r) em.field.E(sm,r),@(t) [0.5 0 t-0.5],[0 1],1)};
for k=1:2
    caught=false;
    try, actions{k}(); catch ex, caught=~isempty(strfind(ex.message,'element 2')); end
    assert(caught,'singularity must identify element 2');
end
bad={@() em.field.V(s,[1 2]),@() em.op.grad(f,[1 2]), ...
    @() em.op.grad(f,p,0),@() em.op.grad(@(r) r,p), ...
    @() em.field.circulation(@(r) r(:,1),seg,[0 1],10), ...
    @() em.field.circulation(FE,@(t) [t t],[0 1],10)};
for k=1:numel(bad)
    caught=false; try, bad{k}(); catch, caught=true; end; assert(caught);
end
% Plot signs/symmetry and arrow direction, independently of perpendicularity.
sd=em.src.merge(em.src.pointCharge(Q,[-0.5 0 0]),em.src.pointCharge(-Q,[0.5 0 0]));
fd=@(r) em.field.V(sd,r); pp=[-1 0.2 0;1 0.2 0;0 1 0]; vv=fd(pp);
assert(vv(1)>0 && vv(2)<0 && abs(vv(1)+vv(2))<1e-13 && abs(vv(3))<1e-13);
[xx,yy]=meshgrid(linspace(-2,2,60)); gridp=[xx(:) yy(:) zeros(numel(xx),1)];
assert(all(isfinite(fd(gridp))));
Edir=em.vec.unit(em.field.E(sd,pp)); assert(all(fd(pp+1e-5*Edir)<fd(pp)));
fig=figure('Visible','off'); ax=axes('Parent',fig);
ch=em.viz.contourV(fd,[-2 2],[-2 2],60,struct('title','Dipole'));
hold on; qh=em.viz.quiver2(@(r) em.field.E(sd,r),[-2 2],[-2 2],16,struct('normalize',true));
assert(isequal(get(ch,'Parent'),ax) && isequal(get(qh,'Parent'),ax));
assert(strcmp(get(get(ax,'XLabel'),'String'),'x (m)') && strcmp(get(get(ax,'YLabel'),'String'),'y (m)'));
assert(~isempty(findall(fig,'Tag','Colorbar')) || ~isempty(findall(fig,'Type','ColorBar')));
em.viz.contourV(fd,[-2 2],[-2 2],60,struct('normalize',true,'logscale',true,'levels',12));
close(fig);
% Default-scale study: report measured optimum, do not assume it.
hs=logspace(-12,-1,23)'; eh=zeros(size(hs)); pe=[1 0.5 0]; ee=FE(pe);
for k=1:numel(hs), eh(k)=norm(-em.op.grad(@(r) em.field.V(s,r),pe,hs(k))-ee)/norm(ee); end
[emin,ii]=min(eh); assert(eh(1)>100*emin && eh(end)>100*emin);
[X,Y]=meshgrid(linspace(0.3,2,40),linspace(-1,1,40)); mp=[X(:) Y(:) zeros(numel(X),1)];
absolute=em.vec.mag(em.field.E(sm,mp)+em.op.grad(f,mp));
tic; largeV=em.field.V(sm,repmat(pe,10000,1)); tv=toc;
tic; largeG=em.op.grad(f,repmat(pe,10000,1)); tg=toc;
assert(isequal(size(largeV),[10000 1]) && isequal(size(largeG),[10000 3]));
assert(tv<2 && tg<2,'10,000-point vectorization timing');
fprintf('  Lab 6: max absolute map error %.3e V/m\n  max closed residual / (max|E| L) %.3e\n',max(absolute),max(ratios(:)));
fprintf('  Square vs exact-tangent edges %.3e V\n  open midpoint errors %s V (order 2)\n',squareDifference,mat2str(er,4));
fprintf('  Fixed line refinement %.3e V\n  growing line potentials %.4f -> %.4f V\n  warning verified\n',abs(diff(vfinite)),vgrow);
fprintf('  Step minimum h=%.3e m, error %.3e\n  default error %.3e\n  10000 points V/grad %.3f/%.3f s\n',hs(ii),emin,eh(hs==1e-5),tv,tg);
fprintf('  test_lab06 checks complete\n');
