function [Data,TrialFrequencies,TrialPhases,TrialDirections,TrialContrasts,timeStepStim,DateTime] = PupilAnalysisToolbox_LoadData(basePath, Protocol,Subject,Date,Session)

% Extract the data path directly from this file.


%*******************
% TO DO -- HANDLE CONTRAST IN THE SAME MANNER AS THE DISTORTION PRODUCT
% The code below was removed from sequential trial analyzer and is placed
% here for now to be integrated
% If the 'Data' struct has a field called 'contrastRelMax', we will
% consider every contrast level as a separate direction. We do that by
% reassigned the direction labels in 'TrialDirections', that way the
% different contrast levels will be handled separately.

% GKA - 8-20-2014
%***************

%dataDirStem = ['/Users/Shared/Matlab/experiments/OneLight/OLPupilDiameter/data/'];

dataDir = fullfile(basePath, Subject, Date, 'MatFiles');
inputFile = [Subject '-' Protocol '-' num2str(Session, '%02.f') '.mat'];

if (exist(fullfile(dataDir, inputFile),'file')==2)
    
    load(fullfile(dataDir, inputFile));
    
    timeStepStim = params.timeStep;
    Data = params.dataStruct;
    TrialFrequencies = [Data.frequencyCarrier];
    
    TrialDirectionCodes = params.theDirections;
    
    % The TrialDirections array is assembled from the modulation cache
    % file names which are more detailed.
    
    TrialDirectionLabels = params.cacheFileName;
    TrialDirections={''};
    for i=1:length(Data)
        TrialDirections(i)=(TrialDirectionLabels(TrialDirectionCodes(i)));
        %TrialDirections(i) = {Data(i).direction};
    end % Assemble a string array of modulation direction cache file names
    
    TrialPhases = [Data.phaseCarrier];
    TrialContrasts = [Data.contrastRelMax];

    % This statemenent detects if the field "modulationMode" exists
    % This detection is needed as some of the older data does not
    % have this field defined
    
    eval('FieldTrue=1;params.dataStruct.modulationMode;','FieldTrue=0;');
    if FieldTrue==1
        
        % Run through the trials. If modulationMode is set and has a value
        % of 'AM", then adjust the data labels for AM mode (distortion product)
        
        CarrierFreqs=[Data.frequencyCarrier];
        EnvelopeFreqs=[Data.frequencyEnvelope];
        EnvelopePhases=[Data.phaseEnvelope];
        
        for i=1:length(Data)
            
            % In distortionProduct studies, the envelope frequency and phase
            % is used to gather and analyze the data. The carrier frequency
            % is effectively a "label" for the direction. The carrier phase
            % is not relevant and is ignored.
            
            % Update the label only if the modulation mode is set to "AM"
            if ~isempty(strfind(TrialDirections{i}, 'Distortion'))
                Data(i).modulationMode = 'AM';
            end
            switch Protocol
                case 'PupillometryDistortionProductTTF'
                    if strcmp(Data(i).modulationMode, 'AM')
                        TrialDirections(i)=strcat(TrialDirections(i),'_',num2str(CarrierFreqs(i), '%03.f'),'Hz');
                        TrialFrequencies(i)=EnvelopeFreqs(i);
                        TrialPhases(i)=EnvelopePhases(i);
                    end
                otherwise
                    if strcmp(Data(i).modulationMode, 'AM')
                        %TrialDirections(i)=strcat(TrialDirections(i),'_',num2str(CarrierFreqs(i), '%03.f'),'Hz');
                        TrialFrequencies(i)=EnvelopeFreqs(i);
                        TrialPhases(i)=EnvelopePhases(i);
                    end
            end
            if strcmp(Data(i).modulationMode, 'pulse')
                TrialPhases(i) = Data(i).phaseRandSec;
            else
                    TrialPhases = TrialPhases/360;
            end
            
            % The following block needs to be used for standard pupil analyses.
            %              if (Data(i).modulationMode == 'AM') & ~isempty(strfind(Protocol, 'TTF'))
            %                 TrialDirections(i)=strcat(TrialDirections(i),'_',num2str(CarrierFreqs(i), '%03.f'),'Hz');
            %                 TrialFrequencies(i)=EnvelopeFreqs(i);
            %                 TrialPhases(i)=EnvelopePhases(i);
            %             end
            %
            
            
        end % Loop through trials looking for AM modulations
    end % Detected the possibility of a modulation mode
    
    % Convert the phase from degree to number of cycles


    
    DateTime = exp.experimentTimeDateString;
    
else % File does not exist; return nulls
    
    Data=[];
    TrialFrequencies=[];
    TrialPhases=[];
    TrialDirections=[];
    timeStepStim=[];
    TrialContrasts = [];
    DateTime = '';
end

end