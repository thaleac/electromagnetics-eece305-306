function s = lineCurrent(I, curve, tspan, N, rule)
%LINECURRENT Discretize a directed filament current.
%   S=EM.SRC.LINECURRENT(I,CURVE,TSPAN,N,RULE) accepts scalar current I
%   in A and the same curve/interval/rules as em.quad.line. Default RULE
%   is 'midpoint'. The source has type 'current', Mx3 positions in m,
%   Mx1 arc-length weights w in m, and Mx3 Idl=I*dr/dt*nodeWeight in A m.
%   Tangents are obtained from em.quad.line, never differentiated again.
%   Stationary quadrature nodes contribute a zero current-length vector.
%
%   Example:
%       c=@(t) [cos(t) sin(t) zeros(size(t))];
%       s=em.src.lineCurrent(1,c,[0 2*pi],400);
%
%   See also EM.QUAD.LINE, EM.SRC.SURFCURRENT, EM.FIELD.H.
if nargin<5||isempty(rule),rule='midpoint';end
if ~isnumeric(I)||~isreal(I)||~isscalar(I)||~isfinite(I)
    error('em:src:lineCurrent:InvalidCurrent','I must be a finite real scalar in A.');
end
[pos,w,tangent]=em.quad.line(curve,tspan,N,rule);
speed=sqrt(sum(tangent.^2,2));factor=zeros(size(w));nz=speed>0;
factor(nz)=I*w(nz)./speed(nz);
s=struct('type','current','pos',pos,'w',w,'Idl',tangent.*factor, ...
    'label','filament current');
end
