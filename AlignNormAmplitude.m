function [xmean ymean xdata ydata] = AlignNormAmplitude(S,varargin)
% ALIGNNORMAMPLITUDE - aligns multiple data sets in y-direction and sets max to zero
%
%   [xmean ymean xdata ydata] = AlignNormAmplitude(S) aligns and normalizes
%   data from struct array S. S should be of the form:
%               S.x:      x-data
%               S.y:      y-data
%               S.weight: weight - logical index indicating whether a
%                         data point should be used for alignment procedure
%                         (i.e. non-significant values can be parsed and
%                         ignored). S.weight should be the same size as S.x
%                         and S.y.
%   For the alignment, y-values with weight = 0 will be set to 0.
%   The output arguments xmean and ymean specify the mean of the original
%   data with the maximum set to 0. Output arguments xdata and ydata are
%   the aligned original data around this mean.
%
%   [xmean ymean xdata ydata] = AlignNormAmplitude(xvalues,yvalues,weights)
%   is similar to the previous call, but uses cell arrays xvalues, yvalues
%   and weights instead of struct array S. Their conversion is as follows:
%
%   xvalues{1}  S(1).x      yvalues{1}  S(1).y      weights{1}  S(1).weight
%	  xvalues{2}  S(2).x      yvalues{2}  S(2).y      weights{2}  S(2).weight
%   ...                     ...                     ...
%   xvalues{n}  S(n).x      yvalues{n}  S(n).y      weights{n}  S(n).weight
%
%   See also AlignNormPhase.
%

% Convert data if necessary
if iscell(S)
    S = cell2struct(S,'x',numel(S));
    [tmp] = cell2struct(varargin{1},'y',numel(S));
    [S.y] = tmp.y;
    [tmp] = cell2struct(varargin{2},'weight',numel(S));
    [S.weight] = tmp.weight;
end

nS = numel(S); % number of datasets

xmin = min([S.x]); % minimum of x
xmax = max([S.x]); % maximum of x
xspac = ceil(log10(xmax-xmin))+3; % increase order of x spacing by 3
narrowx = linspace(xmin,xmax,exp(xspac*log(10))); % interpolation x-values
Xdiff = narrowx(2:end);

WDD = zeros(size(narrowx)); WDD(end) = []; % initialize WDD
Wtot = zeros(size(narrowx)); Wtot(end) = []; % initialize Wtot
for iS = 1:numel(S)
    [Yinterp(iS,:) Ydiff(iS,:) W(iS,:)] = local_diffweight(S(iS),narrowx); % call local_diffweight
    Ydiff(iS,isnan(Ydiff(iS,:))) = 0; % replace Ydiff's NaNs by zeros, otherwise weighted average goes wrong
    Yinterp(iS,isnan(Yinterp(iS,:))) = 0; % replace Yinterp's NaNs by zeros, otherwise weighted average goes wrong
    WDD = WDD+W(iS,:).* Ydiff(iS,:); % weighted diff data
    Wtot = Wtot+W(iS,:); % total weight data
end
warning('off','MATLAB:divideByZero')
MDD = WDD./Wtot; % mean diff data
warning('on','MATLAB:divideByZero')
MDD(isnan(MDD)) = 0; % replace NaNs by zeros

% Output arguments
ymean = cumsum(MDD); % reconstruct average curve
maxymean = max(ymean); % get maximum of curve
ymean = ymean-maxymean; % normalize to zero
xmean = Xdiff;
% Shift and pool original data
xdata = [];
ydata = [];
for iS = 1:numel(S)
    offset(iS) = sum(W(iS,:).*ymean)/sum(W(iS,:))-sum(W(iS,:).*Yinterp(iS,:))/sum(W(iS,:)); % calculate offset
    tmpxdata = S(iS).x;
    tmpydata = S(iS).y; tmpydata(~S(iS).weight) = NaN;
    xdata = [xdata tmpxdata]; % pool xdata
    ydata = [ydata tmpydata+offset(iS)]; % pool ydata
end
ymean(Wtot == 0) = NaN; % set values with zero weight to NaN ('gaps' in data)
[xdata isort] = sort(xdata);
ydata = ydata(isort);

% Set maximum ydata to zero and shift ymean likewise
maxydata = max(ydata);
ydata = ydata-maxydata;
ymean = ymean-maxydata;

ymean = ymean-max(ymean);



function [Yinterp Ydiff W] = local_diffweight(DS,narrowx)
% LOCAL_DIFFWEIGHT - returns interpolated Ydiff and weight

% Set noise floor to zero; valid if applied data is APPLE data
DS.y(~DS.weight) = 0;

% Interpolate y-values
narrowy = interp1(DS.x,DS.y,narrowx);

% Construct weight vector
edges = diff(DS.weight);
on = find(edges > 0)+1; % last element(s) that have weight = 1
off = find(edges < 0); % first element(s) that have weight = 1
if any(on) && any(off)
    if on(1) > off(1) % leading ones in weight
        on = [1 on];
    end
    if on(end) > off(end) % trailing ones in weight
        off = [off numel(DS.weight)];
    end
elseif any(on)
    off = numel(DS.weight);
elseif any(off)
    on = 1;
else
    on = 1;
    off = numel(DS.weight);
end
hdX = diff(DS.x)/2; % mid of difference between DS.x
hdX = [0 hdX 0]; % add first and last indices as extra values; they are only needed for the next two lines
Xon = DS.x(on)-hdX(on); % x-values on
Xoff = DS.x(off)+hdX(off+1); % x-values off
% Actual weights
W = zeros(size(narrowx)); % intialize W (weight)
for ii = 1:numel(Xon) % for all weight-one series
    ind = betwixt(narrowx,Xon(ii),Xoff(ii)); % find indices 
    W(ind) = 1; % set weight to one
end

% Other output args
Yinterp = narrowy(2:end);
Ydiff = diff(narrowy);
W(1) = [];
