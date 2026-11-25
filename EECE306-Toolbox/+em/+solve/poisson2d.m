function [V, info] = poisson2d(G, rho, epsr, mask, Vfix, opts)
%POISSON2D Solve div(eps grad V)=-rho with conservative face fluxes.
%   [V,INFO]=EM.SOLVE.POISSON2D(G,RHO,EPSR,MASK,VFIX,OPTS) returns
%   NyxNx volts. RHO is scalar, NyxNx, or @(X,Y), in C/m^3. EPSR is
%   positive scalar or NyxNx. MASK fixes potentials at true nodes.
%   Unmasked outer faces have zero normal flux (homogeneous Neumann).
%   At least one fixed node is required to remove the constant nullspace.
%
%   OPTS.method: 'direct', 'jacobi', 'gs', or 'sor' (default 'sor');
%   tol=1e-8, maxit=20000, omega=1.9. OPTS.faceMean is 'harmonic'
%   (default) or 'arithmetic' for the Lab 8 comparison. A material jump
%   is represented halfway between adjacent nodes of different materials.
%   INFO contains iter, resHist (maximum potential update in V), time,
%   method, and residual (maximum diagonal-scaled equation residual in V).
%   Direct and iterative methods use the same node-centered flux stencil.
%
%   Example:
%       G=em.solve.grid2d([0 1],[0 1],21,21);
%       m=false(21); m(:,[1 end])=true; v=zeros(21); v(:,1)=1;
%       V=em.solve.poisson2d(G,0,1,m,v,struct('method','direct'));
%
%   See also EM.SOLVE.GRID2D, EM.DERIVE.ENERGYE.
if nargin<6 || isempty(opts), opts=struct(); end
if ~(isstruct(opts)&&isscalar(opts)), error('em:solve:poisson2d:Options','opts must be a scalar struct.'); end
method=option(opts,'method','sor'); tol=option(opts,'tol',1e-8);
maxit=option(opts,'maxit',20000); omega=option(opts,'omega',1.9);
meanType=option(opts,'faceMean','harmonic');
if ~ischar(method)||~ismember(lower(method),{'direct','jacobi','gs','sor'}), error('em:solve:poisson2d:Method','opts.method is invalid.'); end
if ~ischar(meanType)||~ismember(lower(meanType),{'harmonic','arithmetic'}), error('em:solve:poisson2d:Mean','opts.faceMean must be harmonic or arithmetic.'); end
if ~(isnumeric(tol)&&isreal(tol)&&isscalar(tol)&&isfinite(tol)&&tol>0), error('em:solve:poisson2d:Tol','opts.tol must be positive.'); end
if ~(isnumeric(maxit)&&isreal(maxit)&&isscalar(maxit)&&isfinite(maxit)&&maxit>=1&&maxit==floor(maxit)), error('em:solve:poisson2d:Maxit','opts.maxit must be a positive integer.'); end
if ~(isnumeric(omega)&&isreal(omega)&&isscalar(omega)&&isfinite(omega)&&omega>0&&omega<2), error('em:solve:poisson2d:Omega','opts.omega must lie between 0 and 2.'); end
if ~isstruct(G)||~isscalar(G)||~all(isfield(G,{'X','Y','hx','hy'})), error('em:solve:poisson2d:Grid','G must come from grid2d.'); end
sz=size(G.X); ny=sz(1); nx=sz(2);
if numel(sz)~=2||min(sz)<2||~isequal(size(G.Y),sz)||any(~isfinite(G.X(:)))||any(~isfinite(G.Y(:)))||~isscalar(G.hx)||~isscalar(G.hy)||G.hx<=0||G.hy<=0||~isfinite(G.hx)||~isfinite(G.hy)
    error('em:solve:poisson2d:Grid','G must contain finite 2D coordinates and positive spacings.');
end
if isa(rho,'function_handle'), rho=rho(G.X,G.Y); end
R=expand(rho,sz,'rho'); E=expand(epsr,sz,'epsr');
if any(E(:)<=0), error('em:solve:poisson2d:Permittivity','epsr must be positive.'); end
if ~islogical(mask)||~isequal(size(mask),sz)||~any(mask(:)), error('em:solve:poisson2d:Mask','mask must be NyxNx logical with at least one fixed node.'); end
if ~isnumeric(Vfix)||~isreal(Vfix)||~isequal(size(Vfix),sz)||any(~isfinite(Vfix(mask))), error('em:solve:poisson2d:Vfix','Vfix must be NyxNx and finite at mask nodes.'); end
t0=tic;
wx=ones(1,nx); wx([1 end])=0.5; wy=ones(ny,1); wy([1 end])=0.5;
ex=face(E(:,1:end-1),E(:,2:end),meanType).*wy*(G.hy/G.hx);
ey=face(E(1:end-1,:),E(2:end,:),meanType).*wx*(G.hx/G.hy);
ix=reshape(1:ny*nx,sz); p=[reshape(ix(:,1:end-1),[],1);reshape(ix(1:end-1,:),[],1)];
q=[reshape(ix(:,2:end),[],1);reshape(ix(2:end,:),[],1)]; c=[ex(:);ey(:)];
K=sparse([p;q;p;q],[p;q;q;p],[c;c;-c;-c],ny*nx,ny*nx);
b=R.*(wy*wx)*(G.hx*G.hy/em.const.eps0());
u=find(~mask); fixed=find(mask); A=K(u,u); b=b(u)-K(u,fixed)*Vfix(fixed);
V=zeros(sz); V(mask)=Vfix(mask); iter=0; history=[]; residual=0;
if ~isempty(u)
    d=diag(A); x=zeros(numel(u),1);
    if strcmpi(method,'direct')
        x=A\b; iter=1;
    else
        D=spdiags(d,0,numel(u),numel(u)); L=tril(A,-1); U=triu(A,1);
        if strcmpi(method,'sor'), T=D+omega*L; B=omega*U+(omega-1)*D; else, T=D+L; B=U; end
        history=zeros(maxit,1);
        for iter=1:maxit
            old=x;
            if strcmpi(method,'jacobi'), x=x+(b-A*x)./d;
            elseif strcmpi(method,'gs'), x=T\(b-B*x);
            else, x=T\(omega*b-B*x); end
            history(iter)=norm(x-old,inf);
            if history(iter)<tol, break; end
        end
        history=history(1:iter);
        if iter==maxit&&history(end)>=tol, warning('em:solve:poisson2d:NoConvergence','opts.maxit reached before opts.tol.'); end
    end
    V(u)=x; residual=norm((A*x-b)./d,inf);
end
info=struct('iter',iter,'resHist',history,'time',toc(t0),'method',lower(method),'residual',residual);
end
function v=option(o,n,d)
if isfield(o,n), v=o.(n); else, v=d; end
end
function a=expand(a,sz,name)
if isnumeric(a)&&isscalar(a), a=repmat(a,sz); end
if ~isnumeric(a)||~isreal(a)||~isequal(size(a),sz)||any(~isfinite(a(:))), error('em:solve:poisson2d:Input','%s must be finite real scalar or NyxNx.',name); end
end
function f=face(a,b,kind)
if strcmpi(kind,'harmonic'), f=2./(1./a+1./b); else, f=(a+b)/2; end
end
