function videopic(FN)
% videopic - capture and store a picture from video-IN card
%
%   videopic('foo') captures, displays and stores a picture from the
%   video-IN card. There should be no other active connections to
%   the video-IN card, as this will result in an error.
%
%   If videopic is called during an experiment, the picture is stored in 
%   the sub-folder 'pictures' of folder(current(experiment)) and a note
%   with the file's location is added to the experiment's log.
%
%   If there is no active experiment, the user is prompted for a location
%   to store it.
%
%   Example:
%   videopic('test')
%
%   See also imaqhwinfo, videoinput, imaqfind.
%

% Number of frames to use for averaging
nframes = 20; % we want multiple frames to average and remove "the strange black lines"

% Get device information
info = imaqhwinfo('winvideo');
dev_info = imaqhwinfo('winvideo',1);

% Create video input object
vid = videoinput(info.AdaptorName,dev_info.DeviceName,dev_info.SupportedFormats{3}); % dev_info.SupportedFormats: {'RGB24_768x576'  'RGB32_768x576'  'RGB555_768x576'  'UYVY_768x576'  'YUY2_768x576'}

% Set properties of vid
set(vid,'FramesPerTrigger',nframes); % select number of frames to capture

% Run vid if possible
try,
    start(vid);
catch, %#ok<CTCH>
    % Clean up if error was found
    delete(vid);
    return
end

% Retreive data from vid
imgs = getdata(vid);
% Clean up, this is necessary
delete(vid);

% Replace zeros with NaNs
iZs = all(imgs == 0,3); % find all indices where zeros are present
[iZs imgs] = samesize(iZs,imgs); % make iZs same size as imgs
I = double(imgs); % convert image data to double
I(iZs) = NaN; % replace black with NaNs in order to remove "the strange black lines"

% Take average picture
I = nanmean(I,4);
I = uint8(I);

% Show image
imshow(I);

% Save picture to specified location
cexp = current(experiment);
if ~isvoid(cexp),
    dirname = fullfile(folder(cexp),'pictures'); % directory name for pictures
    if ~isdir(dirname)
        mkdir(dirname); % create directory if neccesary in datadir of current experiment
    end
    FFN = fullfilename(FN, dirname, '.tif'); % full file name
    addtolog(cexp,['Picture captured: ' FFN ]);
else,
    % following line is a hack to get (on windows computers) to the user-specific 'My Pictures' dir
    defaultDir = fullfile(fileparts(fileparts(fileparts(fileparts(prefdir)))), 'My Documents', 'My Pictures');
    baseDir = fromsetupFile('LocalSettings', 'defaultPicsDir', '-default', defaultDir);
    [fn, pn] = uiputfile(fullfile(baseDir, ['*.tif']), 'save picture to', fullfilename(FN, baseDir, '.tif'));
    if isequal(0, fn), % user canceled
        FFN = '';
    else,
        FFN = fullfile(pn, fn);
        % remember dir for next call
        tosetupFile('LocalSettings', '-propval', 'defaultPicsDir', pn);
    end
end
if ~isempty(FFN),
    imwrite(I,FFN,'tiff','Compression','none','Resolution',100);
end



