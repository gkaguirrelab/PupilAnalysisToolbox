function [Data,TrialFrequencies,TrialPhases,TrialDirections,AlphaSpacing,DateTime] = DataLoader(Protocol,Subject,Session)

% Extract the data path directly from this file.


%dataDirStem = ['/Users/Shared/Matlab/experiments/OneLight/OLPupilDiameter/data/'];
filePath = fileparts(mfilename('fullpath'));
filePath(strfind(filePath, 'OLPupilDiameter')+length('OLPupilDiameter'):end) = [];
dataDirStem = fullfile(filePath, 'data');

dataDir = fullfile(dataDirStem, Protocol, Subject);
inputFile = [Subject '-' Protocol '-' num2str(Session) '.mat'];

if (exist(fullfile(dataDir, inputFile),'file')==2)
    
    load(fullfile(dataDir, inputFile));
    
    AlphaSpacing = params.alphaSpacing;
    Data = params.dataStruct;
    TrialFrequencies = [0,params.frequencyTrials];  % The initial trial is the adaptation, which has no modulation
    TrialDirectionCodes = [0,params.directionTrials]; % The initial trial is the adaptation, which is no direction
    TrialDirectionLabels = params.cacheFileName;
    TrialDirections={'none'};
    for i=2:length(Data)
        TrialDirections(i)=(TrialDirectionLabels(TrialDirectionCodes(i)));
    end % Assemble a string array of modulation direction cache file names
    
    % Initial operation of the experimental software did not include a
    % variable stimulus phase. The code detects the absence of this
    % parameter and returns a manufactured array of ones, which corresponds
    % to no phase advancement.
    
    if (isfield(params,'phaseTrials'))
        TrialPhases = [1,params.phaseTrials]; % The initial trial is the adaptation, which has no phase delau
    else
        TrialPhases = ones(1,length(Data));
    end
    
    DateTime = exp.experimentTimeDateString;
    
else % File does not exist; return nulls
    
    Data=[];
    TrialFrequencies=[];
    TrialPhases=[];
    TrialDirections=[];
    AlphaSpacing=[];
    DateTime = '';
end

end