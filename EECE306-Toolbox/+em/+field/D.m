function Dv = D(s, r, epsr)
%D Evaluate electric flux density for a discretized charge source.
%   DV = EM.FIELD.D(S, R, EPSR) returns Nx3 values in C/m^2 at the Nx3
%   Cartesian positions R (m). EPSR is a positive scalar relative
%   permittivity, default 1. DV = EPSR*EM.CONST.EPS0()*EM.FIELD.E(S,R).
%
%   Example:
%       s = em.src.pointCharge(1e-9, [0 0 0]);
%       Dv = em.field.D(s, [1 0 0]);
%
%   See also EM.FIELD.E, EM.FIELD.FLUX, EM.CONST.EPS0.

if nargin < 3 || isempty(epsr), epsr = 1; end
if ~(isnumeric(epsr) && isreal(epsr) && isscalar(epsr) ...
        && isfinite(epsr) && epsr > 0)
    error('em:field:D:InvalidPermittivity', ...
        'epsr must be a finite positive real scalar.');
end
Dv = epsr * em.const.eps0() * em.field.E(s, r);
end
