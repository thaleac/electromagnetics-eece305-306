function S = lab10_studies()
%LAB10_STUDIES Evaluate the shared Lab 10 experiments without plotting.
%   Lab10_notebook displays these results, and test_lab10 checks them
%   against the expected values. Run either caller from the toolbox root.
I = 1; a = 0.5; b = 1; c = 1.3; spacing = 0.1;
[gx,gy] = meshgrid((-13:13)*spacing);
% Test membership with integer grid indices to keep the edges symmetric.
ij2 = round(gx/spacing).^2 + round(gy/spacing).^2;
inner = ij2 <= 25;
shield = ij2 >= 100 & ij2 <= 169;
xy = [gx(inner),gy(inner); gx(shield),gy(shield)];
currents = [repmat(I/nnz(inner),nnz(inner),1); repmat(-I/nnz(shield),nnz(shield),1)];
sources = cell(1,size(xy,1));
% Cluster quadrature points near z=0, where all coax measurements lie.
for k = 1:numel(sources)
    sources{k} = straightWire(currents(k),xy(k,:),40,240);
end
coax = em.src.merge(sources{:});
FH = @(r) em.field.H(coax,r);
FB = @(r) em.field.B(coax,r);

rho = [0; 0.12; 0.25; 0.38; 0.47; 0.6; 0.75; 0.9; 1.05; 1.15; 1.25; 1.45; 1.6];
reference = @(r) I*(r<=a).*(r/a).^2 + I*(r>a & r<b) ...
    + I*(r>=b & r<=c).*(c^2-r.^2)/(c^2-b^2);
circulation = zeros(size(rho));
for k = 2:numel(rho)
    circulation(k) = em.field.circulation(FH, circle(rho(k)),[0 2*pi],192);
end
paths = zeros(1,3);
paths(1) = em.field.circulation(FH,circle(0.7),[0 2*pi],96);
corners = 0.65*[-1 -1 0; 1 -1 0; 1 1 0; -1 1 0; -1 -1 0];
for k = 1:4
    start = corners(k,:);
    edge = corners(k+1,:)-start;
    paths(2) = paths(2) + em.field.circulation(FH,@(t) start+t*edge,[0 1],40,'gauss');
end
irregular = @(t) [0.7*(1+0.1*sin(3*t)).*cos(t), ...
    0.68*(1+0.08*cos(2*t)).*sin(t),0*t];
paths(3) = em.field.circulation(FH,irregular,[0 2*pi],128);

segment = straightWire(1,[-1.2 0.15],1.8,180);
ring = em.src.lineCurrent(1,circle(1),[0 2*pi],480);
fluxValues = zeros(2,3);
fluxScale = zeros(2,3);
center = [0.12 -0.09 0.14];
R = 0.55; halfSide = 0.4;
sphere = @(th,ph) center+R*[sin(th).*cos(ph),sin(th).*sin(ph),cos(th)];
disk = @(r,ph) center+[r.*sin(ph),r.*cos(ph),0*r];
faces = {@(u,v)[u v 0*u+halfSide], @(u,v)[v u 0*u-halfSide], ...
    @(u,v)[v 0*u+halfSide u], @(u,v)[u 0*u-halfSide v], ...
    @(u,v)[0*u+halfSide u v], @(u,v)[0*u-halfSide v u]};
magneticSources = {segment,ring};
for k = 1:2
    field = @(r) em.field.B(magneticSources{k},r);
    fluxValues(k,1) = em.field.flux(field,sphere,[0 pi],[0 2*pi],40,72,'gauss');
    curved = em.field.flux(field,sphere,[0 pi/2],[0 2*pi],32,64,'gauss');
    base = em.field.flux(field,disk,[0 R],[0 2*pi],24,64,'gauss');
    fluxValues(k,2) = curved+base;
    for face = 1:6
        position = @(u,v) center+faces{face}(u,v);
        fluxValues(k,3) = fluxValues(k,3) + em.field.flux(field,position, ...
            [-halfSide halfSide],[-halfSide halfSide],24,24,'gauss');
    end
    fluxScale(k,:) = norm(field(center))*[pi*R^2 pi*R^2 (2*halfSide)^2];
end

% The offset sphere intersects the loop without sharing its z=0 symmetry.
cutCenter = [0.92 0.11 0.18]; cutR = 0.46;
cutSphere = @(th,ph) cutCenter+cutR*[sin(th).*cos(ph),sin(th).*sin(ph),cos(th)];
ringField = @(r) em.field.B(ring,r);
cutScale = norm(ringField(cutCenter))*pi*cutR^2;
cutFlux = zeros(1,2);
counts = [160 320];
for k = 1:2
    n = counts(k);
    cutFlux(k) = em.field.flux(ringField,cutSphere,[0 pi],[0 2*pi],n,2*n);
end

J = I/(pi*a^2);
inside = [0 0 0; 0.025 0 0; 0 -0.025 0];
outside = [0.75 0 0; 0 0.76 0; -0.77 0 0];
curlIn = em.op.curl(FH,inside,0.15);
curlOut = em.op.curl(FH,outside,0.015);
divB = em.op.div(FB,[inside;outside],1e-5);
charge = em.src.pointCharge(1e-9,[0.1 0.2 0.3]);
curlE = em.op.curl(@(r) em.field.E(charge,r),[0.45 -0.3 0; -0.75 0.1 0.5; 0.3 0.7 -0.2]);

S = struct('I',I,'spacing',spacing,'edges',[a b c],'counts',[nnz(inner) nnz(shield)], ...
    'rho',rho,'reference',reference,'circulation',circulation,'paths',paths, ...
    'flux',fluxValues,'fluxScale',fluxScale,'cutFlux',cutFlux,'cutScale',cutScale, ...
    'J',J,'curlIn',curlIn,'curlOut',curlOut,'divB',divB,'curlE',curlE);
end

function curve = circle(radius)
curve = @(t) [radius*cos(t),radius*sin(t),0*t];
end

function source = straightWire(I,xy,halfLength,count)
span = asinh(halfLength/0.16);
curve = @(t) [xy(1)+0*t,xy(2)+0*t,0.16*sinh(t)];
source = em.src.lineCurrent(I,curve,[-span span],count);
end
