function data = LoadResultsCSV(theDataDir, theExp, theDataFile, theExpectedHeader)
% Open the file and read the header.  Check that
% what's in the file matches what we expect.
fid = fopen(fullfile(theDataDir,theExp,theDataFile),'r');
tline = fgetl(fid);
theHeader = regexp(tline, '\,', 'split');
for i1 = 1:size(theHeader,2)
    if (~strcmp(theHeader{i1},theExpectedHeader{i1}))
        error('Data file format has changed.  Go check');
    end
end

% Snag the rest of the cells. They all end up as strings
% in a big cell array matrix.
ctr = 0;
while(~feof(fid))
    if ischar(tline)
        ctr = ctr + 1;
        tline = fgetl(fid);
        theCellsAsStrings(ctr,:) = regexp(tline, '\,', 'split');
    else
        break;
    end
end
fclose(fid);

% Pull out the data in a format useful to us.
% Phase angles are in radians [-pi,pi].
% Error is 1 standard deviation.
for i1 = 1:size(theCellsAsStrings,1)
    theDirections{i1} = theCellsAsStrings{i1,1};
    theFrequencies(i1) = str2num(theCellsAsStrings{i1,2});
    theAmplitudes(i1) = str2num(theCellsAsStrings{i1,5});
    thePhasesRadians(i1) = str2num(theCellsAsStrings{i1,6});
    theAmplitudeErr(i1) = str2num(theCellsAsStrings{i1,7});
    thePhaseErrRadians(i1) = str2num(theCellsAsStrings{i1,8});
end

% Sort out data by direction
theMTFNames = unique(theDirections,'stable');
for i1 = 1:length(theMTFNames)
    thisDirection = theMTFNames{i1};
    index = find(ismember(theDirections,thisDirection));
    if (i1 == 1)
        theFrequenciesMTF = theFrequencies(index)';
    else
        theFrequenciesCheck = theFrequencies(index)';
        if (any(theFrequenciesMTF ~= theFrequenciesCheck))
            error('Assumption that all MTFs specified at same frequencies in same order is false');
        end
    end
    theAmplitudeMTF(:,i1) = theAmplitudes(index)';
    thePhaseRadiansMTF(:,i1) = thePhasesRadians(index)';
    theAmplitudeMTFErr(:,i1) = theAmplitudeErr(index)';
    thePhaseMTFErr(:,i1) = thePhaseErrRadians(index)';
end

% Replace the name of the synthetic iso (L/M + S + Mel) with SynIso
theOldName = 'L/M + S + Mel';
theNewName = 'SynIso';
theSynIndex = find(strcmp(theOldName, theMTFNames));
if ~isempty(theSynIndex)
    theMTFNames{theSynIndex} = theNewName;
end

theLog10FrequenciesMTF = log10(theFrequenciesMTF);

% Put everything in a struct
data.theDirections = theDirections;
data.theFrequencies = theFrequencies;
data.theAmplitudes = theAmplitudes;
data.thePhasesRadians = thePhasesRadians;
data.theAmplitudeErr = theAmplitudeErr;
data.thePhaseErrRadians = thePhaseErrRadians;
data.theAmplitudeMTF = theAmplitudeMTF;
data.thePhaseRadiansMTF = thePhaseRadiansMTF;
data.theAmplitudeMTFErr = theAmplitudeMTFErr;
data.thePhaseMTFErr = thePhaseMTFErr;
data.theMTFNames = theMTFNames;
data.theLog10FrequenciesMTF = theLog10FrequenciesMTF;
data.thePath = fullfile(theDataDir,theExp,theDataFile);