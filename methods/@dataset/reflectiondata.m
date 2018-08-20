function [D interpD] = reflectiondata(DS,T)
% dataset/reflectiondata - retrieves laser reflection data for dataset obj.
%
%   D = reflectiondata(DS) returns struct D with fieldnames "Times" and
%   "Samples", containing laser reflection for dataset DS.
%
%   [D interpolatedD] = reflectiondata(DS,T) is similar, but now also returns
%   struct interpolatedD with fieldnames "Times" and "Samples", containing
%   laser reflection for dataset DS interpolated over T.
%
%   See also comport_data/retreivedata
%

% Set default interpolation time
if nargin < 2, T = []; end

% Get record instructions for current dataset
RI = DS.Rec.RecordInstr;

% Get data type index to laser reflection data
iLR = strmatch('Laser reflection', {RI.DataType});

% Error if no laser reflection data available
if isempty(iLR)
    error('No laser reflection data available for specified dataset')
end

% Data fieldname for laser reflection data
dfn = RI(iLR).datafieldname;

% COM-port data
CD = DS.Data.(dfn);

% Retreive and optionally interpolate data
[D interpD] = retrievedata(CD,T);
