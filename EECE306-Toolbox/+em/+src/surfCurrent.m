function s = surfCurrent(K, surface, uspan, vspan, Nu, Nv, rule)
%SURFCURRENT Discretize a vector surface current density.
%   S=EM.SRC.SURFCURRENT(K,SURFACE,USPAN,VSPAN,NU,NV,RULE) uses the
%   area quadrature from em.quad.surf; RULE defaults to 'midpoint'. K is
%   a function handle @(r) returning finite real Mx3 A/m at Mx3 positions,
%   or a 1x3 constant vector. Scalar-only handles are evaluated per row.
%   S.type='current', S.pos is Mx3 m, S.w is Mx1 m^2, and S.Idl=K.*w
%   is Mx3 A m. Quadrature weights are included exactly once in Idl.
%
%   Example:
%       plane=@(u,v) [u v 0*u];
%       s=em.src.surfCurrent(@(r) [1 0 0],plane,[-1 1],[-1 1],20,20);
%
%   See also EM.QUAD.SURF, EM.SRC.LINECURRENT, EM.FIELD.H.
if nargin<7||isempty(rule),rule='midpoint';end
if ~isa(K,'function_handle'),error('em:src:surfCurrent:InvalidK','K must be a function handle returning vectors in A/m.');end
[pos,w]=em.quad.surf(surface,uspan,vspan,Nu,Nv,rule);M=size(pos,1);
try,kv=K(pos);catch,kv=[];end
if ~(isnumeric(kv)&&isreal(kv)&&isequal(size(kv),[M 3]))
    kv=zeros(M,3);
    for k=1:M
        v=K(pos(k,:));
        if ~isnumeric(v)||~isreal(v)||~isequal(size(v),[1 3])
            error('em:src:surfCurrent:InvalidKOutput','K must return a finite real 1x3 or Mx3 vector array.');
        end
        kv(k,:)=v;
    end
end
if any(~isfinite(kv(:))),error('em:src:surfCurrent:InvalidKOutput','K must return finite real vectors in A/m.');end
s=struct('type','current','pos',pos,'w',w,'Idl',kv.*w,'label','surface current');
end
