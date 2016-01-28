% Load the table
T = readtable('/Users/mspits/Dropbox (Aguirre-Brainard Lab)/Melanopsin_Projects/MelLightLevelDependence_Pupillometry/Documents/Datasets.xlsx');
fileList = T(:, 8);
fileList = table2cell(fileList);

% Base folder
baseDir = '/Users/Shared/MATLAB/Experiments/OneLight/OLFlickerSensitivity/data/';

% Target folder
targetDir = '/Users/mspits/Dropbox (Aguirre-Brainard Lab)/Melanopsin_Projects/MelLightLevelDependence_Pupillometry/Data/Archive/SubjectData';

% Make a folder if it doesn't exist
if ~isdir(targetDir)
    mkdir(targetDir);
end

% Iterate over all files
for f = 1:length(fileList)
    if ~isempty(fileList{f})
        RawFile = regexprep(fileList{f}, '/Users/Shared/MATLAB/Experiments/OneLight/OLFlickerSensitivity/data/', '');
        system(['ditto "' fullfile(baseDir, RawFile) '" "' fullfile(targetDir, RawFile) '"']);
    end
end

% Save out the MD5 check sums to control data integrity
currDir = pwd;
cd(targetDir);
system(['find . -type f -exec md5 {} \;>> checksums_archive.md5']);
cd(currDir);

% Expand for cache data

% Expand for calibration data