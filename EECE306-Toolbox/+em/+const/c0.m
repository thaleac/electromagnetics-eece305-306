function c = c0()
%C0 Speed of light in vacuum.
%   C = EM.CONST.C0() returns the speed of light in vacuum in m/s,
%   computed from EM.CONST.MU0 and EM.CONST.EPS0.
%
%   Example:
%       c = em.const.c0()  % approximately 2.998e8 m/s
%
%   See also EM.CONST.EPS0, EM.CONST.MU0.

c = 1 / sqrt(em.const.mu0() * em.const.eps0());
end
