function h = streamlines(F, xBounds, yBounds, seeds, opts)
%STREAMLINES Draw field lines from Cartesian seed points in the z=0 plane.
%   H = EM.VIZ.STREAMLINES(F,XLIM,YLIM,SEEDS,OPTS) returns line handles
%   in seed order. Bounds and Nx2 seed positions are in meters. F takes
%   Nx3 positions and returns finite real Nx3 vectors. Seeds must be in
%   bounds. Choose seeds and bounds away from field singularities.
%   OPTS accepts title, step (m), maxSteps, normalize, and logscale.
%   Positive rescaling of a field leaves the paths unchanged, so the
%   normalize and logscale options do not change the drawn geometry.
%   Tracing stops at a boundary, stagnation, return to the seed, or the
%   step limit. Open paths are traced on both sides of their seeds.
%
%   Example:
%       F = @(r) [-r(:,2), r(:,1), 0*r(:,3)];
%       h = em.viz.streamlines(F,[-1 1],[-1 1],[0.5 0]);
%
%   See also EM.VIZ.QUIVER2.
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isa(F,'function_handle')
    error('em:viz:streamlines:InvalidField','F must be a function handle.');
end
validateattributes(xBounds,{'numeric'},{'real','finite','size',[1 2],'increasing'});
validateattributes(yBounds,{'numeric'},{'real','finite','size',[1 2],'increasing'});
validateattributes(seeds,{'numeric'},{'real','finite','2d','ncols',2});
validateattributes(opts,{'struct'},{'scalar'});
lower = [xBounds(1),yBounds(1)]; upper = [xBounds(2),yBounds(2)];
if any(any(seeds < lower | seeds > upper))
    error('em:viz:streamlines:InvalidSeeds','Seeds must be within the plot bounds.');
end
step = option(opts,'step',min(upper-lower)/240);
maxSteps = option(opts,'maxSteps',3000);
plotTitle = option(opts,'title','Field lines');
validateattributes(step,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(maxSteps,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(option(opts,'normalize',false),{'logical'},{'scalar'});
validateattributes(option(opts,'logscale',false),{'logical'},{'scalar'});

ax = gca;
if ~ishold(ax), cla(ax); end
handles = cell(size(seeds,1),1);
for k = 1:size(seeds,1)
    seed = double(seeds(k,:));
    [xy,closed] = trace(seed,step);
    if ~closed
        reverse = trace(seed,-step);
        xy = [flipud(reverse(2:end,:)); xy];
    end
    colorOrder = get(ax,'ColorOrder');
    color = colorOrder(1+mod(k-1,size(colorOrder,1)),:);
    handles{k} = line(ax,xy(:,1),xy(:,2),'Color',color);
end
h = vertcat(handles{:});
axis(ax,'equal');
xlim(ax,xBounds); ylim(ax,yBounds);
xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); title(ax,plotTitle); grid(ax,'on');

    function [xy,closed] = trace(seed,ds)
        xy = zeros(maxSteps+1,2);
        xy(1,:) = seed;
        closed = false;
        count = 1;
        leftSeed = false;
        for j = 1:maxSteps
            p = xy(count,:);
            k1 = direction(p);
            if ~any(k1), break; end
            k2 = direction(p+ds*k1/2);
            k3 = direction(p+ds*k2/2);
            k4 = direction(p+ds*k3);
            next = p + ds*(k1+2*k2+2*k3+k4)/6;
            if any(next < lower | next > upper) || all(next == p), break; end
            count = count+1;
            xy(count,:) = next;
            distance = norm(next-seed);
            leftSeed = leftSeed || distance > 3*step;
            if leftSeed && distance < step
                xy(count,:) = seed;
                closed = true;
                break;
            end
        end
        xy = xy(1:count,:);
    end

    function v = direction(p)
        if any(p < lower | p > upper)
            v = [0 0];
            return;
        end
        value = F([p 0]);
        if ~(isnumeric(value) && isreal(value) && isequal(size(value),[1 3]) ...
                && all(isfinite(value)))
            error('em:viz:streamlines:InvalidFieldOutput','F must return finite real Nx3 vectors.');
        end
        v = double(value(1:2));
        speed = hypot(v(1),v(2));
        if speed > 0, v = v/speed; end
    end
end

function value = option(opts,name,defaultValue)
value = defaultValue;
if isfield(opts,name), value = opts.(name); end
end
