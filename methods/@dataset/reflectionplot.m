function reflectionplot(DS,varargin)
% dataset/reflectionplot - plots laser reflection data for dataset objects
%
%   reflectionplot(DS,varargin) plots data from dataset DS using LineSpec
%   or other properties specified in varargin.
%
%   See also comport_data/plot.
%

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

% Call plot command
plot(CD,varargin{:})

title(['Laser reflection for dataset ' int2str(DS.ID.iDataset) ' of experiment ' name(DS.ID.Experiment)])