function [p, h] = convergence(errFcn, Nlist, opts)
%CONVERGENCE Measure observed order from errors on successive grids.
%   [P,H] = EM.TEST.CONVERGENCE(ERRFCN, NLIST, OPTS) evaluates ERRFCN(N)
%   for each positive integer in NLIST, fits log(error) against log(N),
%   and returns the negative slope P. H is the log-log plot handle.
%   Set OPTS.plot to false to suppress the plot.
%
%   Example:
%       p = em.test.convergence(@(N) N.^-2, [8 16 32 64]);
%
%   See also POLYFIT, LOGLOG.

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isa(errFcn,'function_handle')
    error('em:test:convergence:InvalidFunction', ...
        'errFcn must be a function handle.');
end
if ~(isnumeric(Nlist) && isreal(Nlist) && isvector(Nlist) ...
        && numel(Nlist) >= 2 && all(isfinite(Nlist(:))) ...
        && all(Nlist(:) > 0) && all(Nlist(:) == floor(Nlist(:))))
    error('em:test:convergence:InvalidNlist', ...
        'Nlist must contain at least two positive integers.');
end
if ~isstruct(opts)
    error('em:test:convergence:InvalidOptions', 'opts must be a struct.');
end
makePlot = true;
if isfield(opts,'plot'), makePlot = opts.plot; end
if ~(islogical(makePlot) && isscalar(makePlot))
    error('em:test:convergence:InvalidPlotOption', ...
        'opts.plot must be a logical scalar.');
end

Nlist = Nlist(:);
err = zeros(size(Nlist));
for k = 1:numel(Nlist)
    value = errFcn(Nlist(k));
    if ~(isnumeric(value) && isreal(value) && isscalar(value) ...
            && isfinite(value) && value > 0)
        error('em:test:convergence:InvalidError', ...
            'errFcn(%d) must return a positive finite scalar.', Nlist(k));
    end
    err(k) = value;
end
fit = polyfit(log(Nlist), log(err), 1);
p = -fit(1);

if makePlot
    h = loglog(Nlist, err, 'o-');
    grid on;
    xlabel('N'); ylabel('error');
    title(sprintf('Observed order p = %.2f', p));
else
    h = [];
end
end
