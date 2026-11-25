function s = merge(varargin)
%MERGE Concatenate source elements for explicit superposition.
%   S = EM.SRC.MERGE(S1, S2, ...) concatenates any number of source
%   structs of the same type. Charge and current sources cannot be merged.
%
%   Example:
%       s = em.src.merge(em.src.pointCharge(1e-9, [-1 0 0]), ...
%                        em.src.pointCharge(1e-9, [ 1 0 0]));
%
%   See also EM.SRC.POINTCHARGE, EM.FIELD.E.

if nargin < 1
    error('em:src:merge:NoSources', ...
        'At least one source argument is required.');
end

first = varargin{1};
validateSource(first, 's1');
sourceType = first.type;
s = first;
s.pos = zeros(0, 3);
s.w = zeros(0, 1);
if strcmp(sourceType, 'charge')
    s.q = zeros(0, 1);
    if isfield(s, 'Idl'), s = rmfield(s, 'Idl'); end
elseif strcmp(sourceType, 'current')
    s.Idl = zeros(0, 3);
    if isfield(s, 'q'), s = rmfield(s, 'q'); end
else
    error('em:src:merge:InvalidType', ...
        's1.type must be ''charge'' or ''current''.');
end
s.label = 'merged source';

for k = 1:nargin
    sk = varargin{k};
    validateSource(sk, sprintf('s%d', k));
    if ~strcmp(sk.type, sourceType)
        error('em:src:merge:TypeMismatch', ...
            'Cannot merge source s%d of type ''%s'' with type ''%s''.', ...
            k, sk.type, sourceType);
    end
    s.pos = [s.pos; sk.pos]; %#ok<AGROW>
    s.w = [s.w; sk.w(:)]; %#ok<AGROW>
    if strcmp(sourceType, 'charge')
        s.q = [s.q; sk.q(:)]; %#ok<AGROW>
    else
        s.Idl = [s.Idl; sk.Idl]; %#ok<AGROW>
    end
end
end

function validateSource(s, name)
if ~isstruct(s) || ~isfield(s, 'type') || ~ischar(s.type)
    error('em:src:merge:InvalidSource', ...
        '%s must be a source struct with a character type field.', name);
end
if ~isfield(s, 'pos') || ~isnumeric(s.pos) || size(s.pos, 2) ~= 3
    error('em:src:merge:InvalidPosition', ...
        '%s.pos must be an Mx3 numeric array.', name);
end
M = size(s.pos, 1);
if ~isfield(s, 'w') || ~isnumeric(s.w) || numel(s.w) ~= M
    error('em:src:merge:InvalidWeights', ...
        '%s.w must contain one numeric weight per source element.', name);
end
if strcmp(s.type, 'charge')
    if ~isfield(s, 'q') || ~isnumeric(s.q) || numel(s.q) ~= M
        error('em:src:merge:InvalidCharge', ...
            '%s.q must contain one numeric charge density per element.', name);
    end
elseif strcmp(s.type, 'current')
    if ~isfield(s, 'Idl') || ~isnumeric(s.Idl) ...
            || ~isequal(size(s.Idl), [M 3])
        error('em:src:merge:InvalidCurrent', ...
            '%s.Idl must be an Mx3 numeric array.', name);
    end
else
    error('em:src:merge:InvalidType', ...
        '%s.type must be ''charge'' or ''current''.', name);
end
end
