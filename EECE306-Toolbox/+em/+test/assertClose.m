function assertClose(actual, expected, tol, msg)
%ASSERTCLOSE  Throw unless ACTUAL and EXPECTED agree to TOL.
%   EM.TEST.ASSERTCLOSE(A, E) compares with a default tolerance of 1e-10.
%   EM.TEST.ASSERTCLOSE(A, E, TOL) sets the tolerance.
%   EM.TEST.ASSERTCLOSE(A, E, TOL, MSG) prepends MSG to the error text.
%
%   The comparison is relative where EXPECTED is nonzero and absolute
%   where it is zero, so that a quantity which should vanish can be
%   checked without dividing by zero.
%
%   Provided with the starter. Do not modify this file.
if nargin < 3 || isempty(tol), tol = 1e-10; end
if nargin < 4, msg = ''; end

a = actual(:);
e = expected(:);
if numel(e) == 1, e = repmat(e, size(a)); end
if numel(a) ~= numel(e)
    error('em:assertClose:sizeMismatch', ...
        '%sactual has %d elements, expected has %d.', tag(msg), numel(a), numel(e));
end

d   = abs(a - e);
rel = e ~= 0;
bad = false(size(a));
bad(rel)  = d(rel) ./ abs(e(rel)) > tol;
bad(~rel) = d(~rel) > tol;

if any(bad)
    k = find(bad, 1);
    if e(k) ~= 0
        detail = sprintf('element %d: got %.12g, expected %.12g, relative error %.3g > %.3g', ...
                         k, a(k), e(k), d(k)/abs(e(k)), tol);
    else
        detail = sprintf('element %d: got %.12g, expected 0, absolute error %.3g > %.3g', ...
                         k, a(k), d(k), tol);
    end
    error('em:assertClose:failed', '%s%s', tag(msg), detail);
end
end

function s = tag(msg)
if isempty(msg), s = ''; else, s = [msg ': ']; end
end
