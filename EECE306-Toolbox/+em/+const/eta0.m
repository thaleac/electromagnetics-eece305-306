function eta = eta0()
%ETA0 - Intrinsic impedance of a vacuum
%   EM.CONST.ETA0() returns the intrinsic impedance in ohms,
%   computed from eps0 and mu0.
%
%   Example:
%       c0 = em.const.c0() % returns approx. 3.00e8
%
%   See also EM.CONST.EPS0, EM.CONST.MU0

    eta = sqrt(em.const.mu0() / em.const.eps0())
end
