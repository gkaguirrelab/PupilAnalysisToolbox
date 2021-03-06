function PupilAnalysisToolbox_MergePupilData(theFileFullPath)
% PupilAnalysisToolbox_MergePupilData
%
% Merges files created on the Windows machine with the experimental file
% generated by OLFlickerSensitivity.
%
% 1/25/16   spitschan   Wrote it.
% 4/6/16    spitschan   Case when there is no input passed.

if ~exist('theFileFullPath', 'var');
   theFileFullPath = []; 
end

if isempty(theFileFullPath)
    %% Open the dialogue to get the file name
    [theFile, thePath] = uigetfile;
    [~, theFileName] = fileparts(theFile);
else
    %% Get the file names from the input
    [thePath, theFileName] = fileparts(theFileFullPath);
    theFile = [theFileName '.mat'];
end

theFilesRaw = dir([thePath filesep theFileName filesep '*.mat']);
nFilesRaw = length(theFilesRaw);

%% Print some input
fprintf('>> Processing %s ...\n', theFile);
fprintf('>> Found %g files ...\n', nFilesRaw);
fprintf('>> File back up in %s...\n', fullfile(thePath, 'unprocessed'));

% Make the back-up directory
if ~exist(fullfile(thePath, 'unprocessed', theFileName), 'dir')
    mkdir(fullfile(thePath, 'unprocessed', theFileName));
end
fprintf('>> Copying files...');
copyfile(fullfile(thePath, theFile), fullfile(thePath, 'unprocessed', theFile));
copyfile(fullfile(thePath, theFileName), fullfile(thePath, 'unprocessed', theFileName));
fprintf('done.\n');

%% Load in the data
load(fullfile(thePath, theFile));

%% Check if the file has already been processed. If yes, return.
if isfield(params, 'mergeState') && params.mergeState
    commandwindow;
    fprintf('*** File has already been merged.\n');
    return;
end

if nFilesRaw ~= length(params.dataStruct)
    fprintf('*** Number of data files does not lign up with what is expected.\n');
    contWithMerging = GetWithDefault(['Continue to merge ' num2str(nFilesRaw) ' raw data files.'], 1);
    if ~contWithMerging
        return;
    end
end


for j = 1:nFilesRaw
    % Load in the raw pupil data
    load(fullfile(thePath, theFileName, theFilesRaw(j).name));
    fprintf('>> Loading file %g / %g ...', j, nFilesRaw);
    dataStructNew(j).diameter = dataStruct.diameter;
    dataStructNew(j).time = dataStruct.time;
    dataStructNew(j).time_inter = dataStruct.time_inter;
    if isempty(dataStruct.time_inter)
        dataStructNew(j).time_inter = 0;
    end
    tmp = allwords(params.cacheFileName{params.theDirections(j)}, '-');
    dataStructNew(j).direction = [tmp{2}];
    tmp0 = params.cacheFileName(params.theDirections);
    if ~isempty(strfind(tmp0{j}, 'Background-45s'));
        phaseSet = [0 90 180 270];
        idx = randi(length(phaseSet), 1);
        dataStructNew(j).frequencyCarrier = 0.1;
        dataStructNew(j).phaseCarrier = phaseSet(idx);
        dataStructNew(j).modulationMode = 'FM';
    end
    dataStructNew(j).modulationMode = 'pulse';
    if strcmp(dataStructNew(j).modulationMode, 'pulse')
        theShifts = 0:5;
        dataStructNew(j).frequencyCarrier = -1;
        dataStructNew(j).frequencyEnvelope = -1;
        dataStructNew(j).phaseCarrier = -1;
        dataStructNew(j).phaseEnvelope = -1;
        dataStructNew(j).phaseRandSec = theShifts(params.thePhaseIndices(j));
    end
    % The following section is commented out; this will need to be
    % addressed at a later point.
    %     if ~isempty(strfind(protocol, 'CRF'))
    %         dataStructNew(j).direction = [dataStructNew(j).direction '_CRF_' num2str(abs(dataStructNew(j).contrastRelMax)*42, '%02.fpct')];
    %     elseif ~isempty(strfind(protocol, 'DoublePulse'))
    %         if dataStructNew(j).contrastRelMax >= 0
    %             dataStructNew(j).direction = [dataStructNew(j).direction '_DoublePulse_' num2str(abs(dataStructNew(j).contrastRelMax)*42, '%02.fpct')];
    %         elseif dataStructNew(j).contrastRelMax <0
    %             dataStructNew(j).direction = [dataStructNew(j).direction '_DoublePulse_-' num2str(abs(dataStructNew(j).contrastRelMax)*42, '%02.fpct')];
    %         end
    %     else
    dataStructNew(j).contrastRelMax = 1;
    %end
    
    % Save out some of the raw pupil stuff
    dataStructNew(j).rawTimeStamps = pupilData.timeStamps;
    dataStructNew(j).rawPupilDiameter = pupilData.pupilDiameter;
    dataStructNew(j).rawMmPositions = pupilData.mmPositions;
    dataStructNew(j).rawFickPositions = pupilData.fickPositions;
    dataStructNew(j).rawHelmholtzPositions = pupilData.helmholtzPositions;
    fprintf('done.\n');
end

params.dataStruct = dataStructNew;
params.mergeState = true;
save(fullfile(thePath, theFile), 'params', 'exp');
fprintf('*** All merging finished.\n');