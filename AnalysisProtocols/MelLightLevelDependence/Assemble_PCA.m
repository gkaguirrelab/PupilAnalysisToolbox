basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sPulse_';

M_LMS = [];
M_LMS_avg = [];
M_Mel = [];

%% Load CSV files
ndVal = 'ND20';
theSubjects = {'G092815Ax' 'J092915Rx' 'M092515Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+_TimeSeries.csv']));
    M_LMS = [M_LMS tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+_TimeSeries.csv']));
    M_Mel = [M_Mel tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS_avg = [M_LMS_avg tmp(:, 2)];
end


%% Load CSV files
ndVal = 'ND15';
theSubjects = {'G100515Ax' 'J100615Rx' 'M100615Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+_TimeSeries.csv']));
    M_LMS = [M_LMS tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+_TimeSeries.csv']));
    M_Mel = [M_Mel tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS_avg = [M_LMS_avg tmp(:, 2)];
end

%% Load CSV files
ndVal = 'ND10';
theSubjects = {'G092815Ax' 'J092915Rx' 'M092515Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+_TimeSeries.csv']));
    M_LMS = [M_LMS tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+_TimeSeries.csv']));
    M_Mel = [M_Mel tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS_avg = [M_LMS_avg tmp(:, 2)];
end




%% Load CSV files
ndVal = 'ND05';
theSubjects = {'G100515Ax' 'J100615Rx' 'M100615Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+_TimeSeries.csv']));
    M_LMS = [M_LMS tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+_TimeSeries.csv']));
    M_Mel = [M_Mel tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS_avg = [M_LMS_avg tmp(:, 2)];
end



%% Load CSV files
ndVal = 'ND00';
theSubjects = {'G092815Ax' 'J092915Rx' 'M092515Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+_TimeSeries.csv']));
    M_LMS = [M_LMS tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+_TimeSeries.csv']));
    M_Mel = [M_Mel tmp];
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS_avg = [M_LMS_avg tmp(:, 2)];
end

%%

meanLMS = mean(M_LMS_avg, 2);

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU]   = pca([M_Mel M_LMS]');
nComponents = 10;
subplot(1, 2, 1);
plot(1:nComponents, EXPLAINED(1:nComponents), '-ok', 'MarkerFaceColor', 'k');
xlabel('Components');
ylabel('Variance explained [%]');
pbaspect([1 1 1]);

subplot(1, 2, 2);
plot(MU, '-k', 'LineWidth', 1.5); hold on;
plot(COEFF(:, 1)); hold on;

plot(COEFF(:, 3))

plot(COEFF(:, 2))
pbaspect([1 1 1]);