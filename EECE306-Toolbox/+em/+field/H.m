function Hv = H(s, r)
%H Biot-Savart field of a finite discretized current source.
%   HV=EM.FIELD.H(S,R) returns Nx3 A/m at finite real Nx3 positions R
%   in meters. S.type is 'current'; S.pos is Mx3 and S.Idl is Mx3
%   directed current-length (A m), already including quadrature weights.
%   HV=sum(Idl cross (r-pos)/|r-pos|^3)/(4*pi). No extra S.w factor
%   is applied. Observation/source pairs are vectorized in bounded chunks.
%   Positions within 1e-12 times the coordinate scale of an element raise
%   a singularity error naming the observation and source element.
%
%   Example:
%       c=@(t) [cos(t) sin(t) zeros(size(t))];
%       s=em.src.lineCurrent(1,c,[0 2*pi],400);
%       H=em.field.H(s,[0 0 0]); % approximately [0 0 0.5] A/m
%
%   See also EM.FIELD.E, EM.FIELD.B, EM.SRC.LINECURRENT.
if ~isstruct(s)||~isscalar(s)||~isfield(s,'type')||~strcmp(s.type,'current')
    error('em:field:H:InvalidSource','s must be a scalar current source struct.');
end
if ~isfield(s,'pos')||~isnumeric(s.pos)||~isreal(s.pos)||~ismatrix(s.pos)||size(s.pos,2)~=3||any(~isfinite(s.pos(:)))
    error('em:field:H:InvalidSource','s.pos must be finite real Mx3.');
end
M=size(s.pos,1);
if ~isfield(s,'Idl')||~isnumeric(s.Idl)||~isreal(s.Idl)||~isequal(size(s.Idl),[M 3])||any(~isfinite(s.Idl(:)))
    error('em:field:H:InvalidSource','s.Idl must be finite real Mx3.');
end
if ~isfield(s,'w')||~isnumeric(s.w)||~isreal(s.w)||numel(s.w)~=M||any(~isfinite(s.w(:)))
    error('em:field:H:InvalidSource','s.w must contain M finite real weights.');
end
if ~isnumeric(r)||~isreal(r)||~ismatrix(r)||size(r,2)~=3||any(~isfinite(r(:)))
    error('em:field:H:InvalidPosition','r must be a finite real Nx3 position array.');
end
Hv=zeros(size(r)); if isempty(r)||M==0,return;end
chunk=max(1,floor(1e6/M));
tol=1e-12*max([1;abs(r(:));abs(s.pos(:))]);
j=reshape(s.Idl,[1 M 3]);
for first=1:chunk:size(r,1)
    last=min(first+chunk-1,size(r,1)); nc=last-first+1;
    d=reshape(r(first:last,:),[nc 1 3])-reshape(s.pos,[1 M 3]);
    R2=sum(d.^2,3);[obs,elem]=find(R2<=tol^2,1);
    if ~isempty(elem)
        error('em:field:H:Singularity','Observation point %d lies at or too near source element %d; H is singular.',first+obs-1,elem);
    end
    % Spell out Idl cross d to avoid expanding the source array in memory.
    crossJ=cat(3,j(:,:,2).*d(:,:,3)-j(:,:,3).*d(:,:,2), ...
        j(:,:,3).*d(:,:,1)-j(:,:,1).*d(:,:,3), ...
        j(:,:,1).*d(:,:,2)-j(:,:,2).*d(:,:,1));
    Hv(first:last,:)=reshape(sum(crossJ./(R2.*sqrt(R2)),2),[nc 3]);
end
Hv=Hv/(4*pi);
end
