function rho = mms2d(Vexact, epsr)
%MMS2D Build the charge density for a manufactured 2D potential.
%   RHO = EM.TEST.MMS2D(VEXACT,EPSR) returns a function handle @(X,Y)
%   for rho = -eps0*epsr*Laplacian(Vexact). VEXACT must accept equally
%   sized X and Y arrays. EPSR is the positive uniform relative
%   permittivity used by the Lab 7 Poisson stencil.
%
%   Example:
%       V = @(X,Y) sin(pi*X).*sin(pi*Y);
%       rho = em.test.mms2d(V, 1);
%
%   See also EM.SOLVE.POISSON2D, EM.OP.LAPLACIAN.

if ~isa(Vexact, 'function_handle')
    error('em:test:mms2d:InvalidPotential', ...
        'Vexact must be a function handle @(X,Y).');
end
if ~(isnumeric(epsr) && isreal(epsr) && isscalar(epsr) ...
        && isfinite(epsr) && epsr > 0)
    error('em:test:mms2d:InvalidPermittivity', ...
        'epsr must be a finite positive scalar.');
end

e0 = em.const.eps0();
h = eps^(1/4);
rho = @(X,Y) source(Vexact, e0*epsr, X, Y, h);
end

function R = source(Vexact, epsValue, X, Y, h)
if ~(isnumeric(X) && isreal(X) && isnumeric(Y) && isreal(Y) ...
        && isequal(size(X), size(Y)) && all(isfinite(X(:))) ...
        && all(isfinite(Y(:))))
    error('em:test:mms2d:InvalidCoordinates', ...
        'X and Y must be equally sized finite real arrays.');
end
V0 = Vexact(X,Y);
Vxp = Vexact(X+h,Y); Vxm = Vexact(X-h,Y);
Vyp = Vexact(X,Y+h); Vym = Vexact(X,Y-h);
values = {V0,Vxp,Vxm,Vyp,Vym};
for k = 1:numel(values)
    value = values{k};
    if ~(isnumeric(value) && isreal(value) && isequal(size(value),size(X)) ...
            && all(isfinite(value(:))))
        error('em:test:mms2d:InvalidPotentialOutput', ...
            'Vexact must return a finite real array matching X and Y.');
    end
end
R = -epsValue * ((Vxp - 2*V0 + Vxm) + (Vyp - 2*V0 + Vym)) / h^2;
end
