function Bv = B(s, r, mur)
%B Magnetic flux density from a discretized current source.
%   BV=EM.FIELD.B(S,R,MUR) returns Nx3 tesla at Nx3 positions R (m),
%   as mur*em.const.mu0()*em.field.H(S,R). MUR is a positive finite
%   scalar relative permeability, default 1. This is a uniform-medium
%   conversion, not a solver for spatial permeability interfaces.
%
%   Example:
%       s=em.src.lineCurrent(1,@(t) [cos(t) sin(t) 0*t],[0 2*pi],400);
%       B=em.field.B(s,[0 0 0]);
%
%   See also EM.FIELD.H, EM.CONST.MU0, EM.FIELD.D.
if nargin<3||isempty(mur),mur=1;end
if ~isnumeric(mur)||~isreal(mur)||~isscalar(mur)||~isfinite(mur)||mur<=0
    error('em:field:B:InvalidMur','mur must be a positive finite real scalar.');
end
Bv=mur*em.const.mu0()*em.field.H(s,r);
end
