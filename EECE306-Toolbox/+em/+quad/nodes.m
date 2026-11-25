function [t, w] = nodes(a, b, N, rule)
%NODES Quadrature nodes and weights on a finite interval.
%   [T, W] = EM.QUAD.NODES(A, B, N, RULE) returns Nx1 nodes and Nx1
%   weights on [A,B]. RULE is 'midpoint', 'trapz', 'simpson', or 'gauss'.
%   The default rule is 'midpoint'.
%
%   For an even number of Simpson nodes, the final three subintervals use
%   Simpson's 3/8 rule and the preceding even number of subintervals use
%   Simpson's 1/3 rule. This preserves N nodes and fourth-order accuracy.
%
%   Example:
%       [t,w] = em.quad.nodes(0, 1, 5, 'simpson');
%       integralEstimate = sum(t.^3 .* w);
%
%   See also EM.QUAD.LINE, EM.QUAD.SURF, EM.QUAD.VOL.

if nargin < 4 || isempty(rule), rule = 'midpoint'; end
if ~(isnumeric(a) && isreal(a) && isscalar(a) && isfinite(a))
    error('em:quad:nodes:InvalidA', 'a must be a finite real scalar.');
end
if ~(isnumeric(b) && isreal(b) && isscalar(b) && isfinite(b) && b > a)
    error('em:quad:nodes:InvalidB', ...
        'b must be a finite real scalar greater than a.');
end
if ~(isnumeric(N) && isreal(N) && isscalar(N) && isfinite(N) ...
        && N == floor(N) && N >= 1)
    error('em:quad:nodes:InvalidN', 'N must be a positive integer scalar.');
end
if isstring(rule) && isscalar(rule), rule = char(rule); end
if ~ischar(rule)
    error('em:quad:nodes:InvalidRule', 'rule must be character text.');
end

rule = lower(strtrim(rule));
switch rule
    case 'midpoint'
        h = (b-a)/N;
        t = a + ((0:N-1)' + 0.5)*h;
        w = repmat(h, N, 1);

    case 'trapz'
        if N < 2
            error('em:quad:nodes:TrapzNodeCount', ...
                'N must be at least 2 for the trapezoid rule.');
        end
        h = (b-a)/(N-1);
        t = linspace(a, b, N)';
        w = repmat(h, N, 1);
        w([1 end]) = h/2;

    case 'simpson'
        if N < 3
            error('em:quad:nodes:SimpsonNodeCount', ...
                'N must be at least 3 for the Simpson rule.');
        end
        h = (b-a)/(N-1);
        t = linspace(a, b, N)';
        w = zeros(N,1);
        if mod(N,2) == 1
            w(1:2:end) = 2;
            w(2:2:end-1) = 4;
            w([1 end]) = 1;
            w = (h/3)*w;
        else
            if N > 4
                last13 = N-3;
                w(1:2:last13) = w(1:2:last13) + 2;
                w(2:2:last13-1) = w(2:2:last13-1) + 4;
                w(1) = w(1) - 1;
                w(last13) = w(last13) - 1;
                w = (h/3)*w;
            end
            w(N-3:N) = w(N-3:N) + (3*h/8)*[1; 3; 3; 1];
        end

    case 'gauss'
        if N == 1
            x = 0;
            wg = 2;
        else
            k = (1:N-1)';
            beta = k ./ sqrt(4*k.^2 - 1);
            J = diag(beta,1) + diag(beta,-1);
            [V,D] = eig(J);
            [x,order] = sort(diag(D));
            V = V(:,order);
            wg = 2*(V(1,:)'.^2);
        end
        t = (a+b)/2 + ((b-a)/2)*x;
        w = ((b-a)/2)*wg;

    otherwise
        error('em:quad:nodes:UnknownRule', ...
            'rule must be ''midpoint'', ''trapz'', ''simpson'', or ''gauss''.');
end

t = t(:);
w = w(:);
end
