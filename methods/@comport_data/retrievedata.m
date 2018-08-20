function [D interpD] = retreivedata(CD,T)
% comport_data/retrievedata - retrieves data for comport_data objects
%
%   D = retrievedata(CD) returns struct D with fieldnames "Times" and
%   "Samples", containing data retrieved from comport_data object CD.
%
%   [D interpolatedD] = retrievedata(CD,T) is similar, but now also returns
%   struct interpolatedD with fieldnames "Times" and "Samples", containing
%   data from comport_data object CD interpolated over T.
%

% Set default interpolation time
if nargin < 2, T = []; end

% Access actual data
D.Times = CD.Data.Times;
D.Samples = CD.Data.Samples;

% Interpolate data if asked
interpD.Times = T;
interpD.Samples = nan(size(T));
if ~isempty(T)
    interpD.Samples = interp1(D.Times,D.Samples,T);
end
