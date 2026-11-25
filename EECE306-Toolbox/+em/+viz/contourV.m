function h = contourV(f, xBounds, yBounds, N, opts)
%CONTOURV Filled scalar-field contours on the z=0 plane in current axes.
%   H = EM.VIZ.CONTOURV(F,XLIM,YLIM,N,OPTS) samples F at an NxN grid
%   and returns the contour handle. F accepts Nx3 and returns Nx1.
%   OPTS defaults to struct(); title defaults to 'Electric potential'.
%   normalize (logical, default false) divides by max(abs(V)); an all-zero
%   field stays zero. logscale (logical, default false) applies signed
%   log10(1+abs(V)), after normalization if enabled, retaining dipole signs.
%   These monotone transforms preserve level sets. The colorbar states the
%   transform; raw values use volts. levels defaults to 20 (positive
%   integer count or finite increasing vector). Existing axes support
%   hold on and EM.VIZ.QUIVER2 overlays; no new figure is created.
%
%   Example:
%       s=em.src.pointCharge(1e-9,[0 0 0]);
%       em.viz.contourV(@(r) em.field.V(s,r),[-2 2],[-2 2],60,struct());
%
%   See also EM.FIELD.V, EM.VIZ.QUIVER2.
if nargin<5 || isempty(opts), opts=struct(); end
if ~isa(f,'function_handle'), error('em:viz:contourV:InvalidField','f must be a function handle.'); end
for b={xBounds,yBounds}
    a=b{1};
    if ~(isnumeric(a) && isreal(a) && isequal(size(a),[1 2]) && all(isfinite(a)) && a(2)>a(1))
        error('em:viz:contourV:InvalidBounds','xlim and ylim must be finite increasing 1x2 arrays.');
    end
end
if ~(isnumeric(N) && isreal(N) && isscalar(N) && isfinite(N) && N>=2 && N==floor(N))
    error('em:viz:contourV:InvalidN','N must be an integer scalar >=2.');
end
if ~(isstruct(opts) && isscalar(opts)), error('em:viz:contourV:InvalidOptions','opts must be a scalar struct.'); end
normalize=option(opts,'normalize',false); logscale=option(opts,'logscale',false);
plotTitle=option(opts,'title','Electric potential'); levels=option(opts,'levels',20);
if ~(islogical(normalize) && isscalar(normalize) && islogical(logscale) && isscalar(logscale))
    error('em:viz:contourV:InvalidOptions','opts.normalize and opts.logscale must be logical scalars.');
end
if ~(ischar(plotTitle) || (isstring(plotTitle) && isscalar(plotTitle)))
    error('em:viz:contourV:InvalidTitle','opts.title must be character text.');
end
if ~(isnumeric(levels) && isreal(levels) && isvector(levels) && all(isfinite(levels)) && ...
        ((isscalar(levels) && levels>=1 && levels==floor(levels)) || (numel(levels)>1 && all(diff(levels)>0))))
    error('em:viz:contourV:InvalidLevels','opts.levels must be a positive count or increasing vector.');
end
[X,Y]=meshgrid(linspace(xBounds(1),xBounds(2),N),linspace(yBounds(1),yBounds(2),N));
values=f([X(:) Y(:) zeros(numel(X),1)]);
if ~(isnumeric(values) && isreal(values) && isequal(size(values),[numel(X) 1]) && all(isfinite(values(:))))
    error('em:viz:contourV:InvalidFieldOutput','f must return finite real Nx1 values.');
end
label='V (V)';
if normalize
    scale=max(abs(values)); if scale>0, values=values/scale; end
    label='V / max|V|';
end
if logscale, values=sign(values).*log10(1+abs(values)); label=['signed log10(1+abs(' label '))']; end
[~,h]=contourf(X,Y,reshape(values,size(X)),levels);
axis equal; xlim(xBounds); ylim(yBounds); xlabel('x (m)'); ylabel('y (m)');
title(plotTitle); cb=colorbar; ylabel(cb,label); grid on;
end
function v=option(opts,name,default)
if isfield(opts,name), v=opts.(name); else, v=default; end
end
