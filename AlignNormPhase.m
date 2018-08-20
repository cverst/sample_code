function [xmean ymean xdata ydata] = AlignNormPhase(S,varargin)
% ALIGNNORMPHASE - aligns multiple phase data sets and sets -1<origin<=0
%
%   [xmean ymean xdata ydata] = AlignNormPhase(S) aligns and normalizes
%   phase data from struct array S. S should be of the form:
%               S.x:      x-data
%               S.y:      y-data
%               S.weight: weight - logical index indicating whether a
%                         data point should be used for alignment procedure
%                         (i.e. non-significant values can be parsed and
%                         ignored). S.weight should be the same size as S.x
%                         and S.y.
%   For the alignment, y-values with weight = 0 will be ignored.
%   The output arguments xmean and ymean specify the mean of the original
%   data with the origin set between -1 and 0 cylces. Output arguments
%   xdata and ydata are the original data points around this mean.
%
%   [xmean ymean xdata ydata] = AlignNormPhase(xvalues,yvalues,weights) is
%   similar to the previous call, but uses cell arrays xvalues, yvalues and
%   weights instead of struct array S. Their conversion is as follows:
%
%   xvalues{1}  S(1).x      yvalues{1}  S(1).y      weights{1}  S(1).weight
%	  xvalues{2}  S(2).x      yvalues{2}  S(2).y      weights{2}  S(2).weight
%   ...                     ...                     ...
%   xvalues{n}  S(n).x      yvalues{n}  S(n).y      weights{n}  S(n).weight
%
%   See also AlignNormAmplitude.
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

for iS = 1:numel(S)
    [Yinterp(iS,:) W(iS,:)] = local_interpYweight(S(iS),narrowx); % call local_interpYweight
    Wnan = ones(size(W(iS,:))); Wnan(~W(iS,:)) = NaN;
    WD(iS,:) = Wnan.*exp(i*2*pi*Yinterp(iS,:)); % weighted data
end
MD = nanmean(WD,1); % weighted mean data

% Output arguments
[junk nanindx] = denan(MD);
dymean = cunwrap(cangle(MD(nanindx))); % set phase to cycles and unwrap
ymean = nan(size(MD));
ymean(nanindx) = dymean;
if numel(dymean) ~= 0
    startymean = ceil(dymean(1)); % get start y-value of curve
else
    startymean = NaN;
end
ymean = ymean-startymean; % set start between -1 and 0
xmean = narrowx(2:end);

% Shift and pool original data
xdata = [];
ydata = [];
for iS = 1:numel(S)
    xsign = S(iS).x(S(iS).weight); % significant elements of x
    ysign = S(iS).y(S(iS).weight); % significant elements of y
    if ~isempty(xsign)
        warning('off','MATLAB:interp1:NaNinY')
        yi = interp1(xmean,ymean,xsign);
        ikeep = find(~isnan(ymean));
        if isnan(yi(1))
            yi(1) = ymean(ikeep(1));
        end
        if isnan(yi(end))
            yi(end) = ymean(ikeep(end));
        end
        tempoffset = -round(ysign-yi);
        warning('on','MATLAB:interp1:NaNinY')
        offset = zeros(size(S(iS).weight));
        offset(find(S(iS).weight == true)) = tempoffset;
    else
        offset = NaN;
    end
    tmpxdata = S(iS).x;
    tmpydata = S(iS).y; tmpydata(~S(iS).weight) = NaN;
    xdata = [xdata tmpxdata]; % pool xdata
    ydata = [ydata tmpydata+offset]; % pool ydata
end
[xdata isort] = sort(xdata);
ydata = ydata(isort);


function [Yinterp W] = local_interpYweight(DS,narrowx)
% LOCAL_INTERPYWEIGHT - returns interpolated Y and weight

% Remove non-significant data for interpolation
Y = DS.y;
X = DS.x;
Y(~DS.weight) = [];
X(~DS.weight) = [];

% Interpolate y-values without non-significant data
if isempty(X)
    narrowy = nan(size(narrowx));
elseif numel(X) > 1
    narrowy = interp1(X,Y,narrowx);
else
    narrowy = ones(size(narrowx))*Y;
end

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
hdX = diff(DS.x)/2; % mid of difference between DS.x (half diff of X)
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
W(1) = [];
