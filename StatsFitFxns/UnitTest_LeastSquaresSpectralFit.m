clc; clear all; close all;

theUnitTest = figure;
% Cosine wave has a phase of 90
sampling_frequency = 20;
t = 0:1/sampling_frequency:100;
frequency = 0.1;
data = -1*cos(2*pi*frequency*t);

[~,FitPhaseNew,~] = LeastSquaresSpectralFit(data, frequency, sampling_frequency);
[~,FitPhaseStutter,~] = LeastSquaresSpectralFit_StutterError(data, frequency, sampling_frequency);

theSignalPhase = -90;
disp(['Cosine wave']);
disp(['- Signal phase: ' num2str(wrapTo180(theSignalPhase))]);
disp(['- Calculated phase (new): ' num2str(FitPhaseNew)]);
disp(['- Calculated phase (new, stutter-corrected): ' num2str(FitPhaseStutter)]);

% Sine wave has a phase of 0
sampling_frequency = 20;
t = 0:1/sampling_frequency:100;
frequency = 0.1;
data = -1*sin(2*pi*frequency*t);

[~,FitPhaseNew,BestFit] = LeastSquaresSpectralFit(data, frequency, sampling_frequency);
[~,FitPhaseStutter,BestFitStutter] = LeastSquaresSpectralFit_StutterError(data, frequency, sampling_frequency);
theSignalPhase = 0;
disp(['Sine wave']);
disp(['- Signal phase: ' num2str(wrapTo180(theSignalPhase))]);
disp(['- Calculated phase (new): ' num2str(wrapTo180((FitPhaseNew)))]);
disp(['- Calculated phase (new, stutter-corrected): ' num2str(wrapTo180((FitPhaseStutter)))]);

% The phase
thePhases = [-pi:0.1:pi];

sampling_frequency = 20;
t = 0:1/sampling_frequency:100;
frequency = 0.1;
figure;
for i = 1:length(thePhases)
    data = -1*sin(2*pi*frequency*t + thePhases(i));
    
    [~,FitPhaseNew(i)] = LeastSquaresSpectralFit(data, frequency, sampling_frequency);
    [~,FitPhaseStutter(i)] = LeastSquaresSpectralFit_StutterError(data, frequency, sampling_frequency);
    
    plot(i, thePhases(i), 'or'); hold on;
    plot(i, FitPhaseNew(i), 'xk');
    
end

figure;
plot(-1*sin(2*pi*frequency*t + thePhases(28)), '-b'); hold on;
plot(-1*sin(2*pi*frequency*t), '--k');
legend('Delay');

% Sine wave with the a fraction of signal cycle chopped
sampling_frequency = 20;
t = 0:1/sampling_frequency:50;
frequency = .1;
data = -1*sin(2*pi*frequency*t);
choppedCycle = 0.6;
choppedTime = choppedCycle* (1/frequency);
choppedSample = round(choppedTime * sampling_frequency);
data(1:choppedSample) = NaN;
[~,FitPhaseNew,BestFit] = LeastSquaresSpectralFit(data, frequency, sampling_frequency);
[~,FitPhaseStutter,BestFitStutter] = LeastSquaresSpectralFit_StutterError(data, frequency, sampling_frequency);
theSignalPhase = 0;
disp(['Sine wave']);
disp(['- Signal phase: ' num2str(wrapTo180(theSignalPhase))]);
disp(['- Calculated phase (new): ' num2str(wrapTo180((FitPhaseNew)))]);
disp(['- Calculated phase (new, stutter-corrected): ' num2str(wrapTo180((FitPhaseStutter)))]);


% Phase shift from 0 to 1 cycle
originalSignal = -1*sin(2*pi*frequency*t);
phaseShiftInSecond = [0 0.25 0.5 0.75 1]*(1/frequency);
for i = 1:length(phaseShiftInSecond);
    thePhase = phaseShiftInSecond(i);
    signalShift = -1*sin(2*pi*frequency*t - phaseShiftInSecond(i)*2*pi);
    signalShift(1:choppedSample) = NaN;
    
    [~,FitPhaseNew,BestFit] = LeastSquaresSpectralFit(signalShift, frequency, sampling_frequency);
    [~,FitPhaseStutter,~] = LeastSquaresSpectralFit_StutterError(signalShift, frequency, sampling_frequency);
    theSignalPhase(i) = thePhase;
    theFitPhaseNew(i) = FitPhaseNew;
    theFitPhaseStutter(i) = FitPhaseStutter;
    
    figure;
    plot(t,originalSignal, 'Color', [0.5 0.5 0.5]); hold on;
    plot(t,signalShift, 'r');
end

% Convert the phase to number of cycles
theSignalPhase = theSignalPhase * frequency;
theFitPhaseNew = unwrap(theFitPhaseNew);
theFitPhaseNew = theFitPhaseNew/(2*pi);

theUnitTest = figure;
plot(theSignalPhase, 'xk','MarkerSize',10); hold on
plot(theFitPhaseNew, 'ob','MarkerSize',10); hold on
%plot(wrapTo2Pi(-theFitPhaseNew), '-r');
legend('Simulated phase', 'Fit phase', 'Fit phase *-1');
ylabel('Phase [cycle]');
title('Phase');

saveFileName = 'UnitTest_LeastSquaresSpectralFit';
set(theUnitTest, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(theUnitTest, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(theUnitTest, saveFileName, 'pdf');

%% Test the normalization error
a = rand(1,3);
b = rand(1,3);

% Right theta
A = sqrt(a.^2 + b.^2);
theta = atan2(b,a);
complexSum = sum(A.*exp(1i*theta));
thetaSum = phase(complexSum);
thetaNorm = theta - thetaSum;

% Wrong theta
Ahat = sqrt(a.^2 + b.^2);
thetaHat = atan2(a,b);
complexSumHat = sum(Ahat.*exp(1i*thetaHat));
thetaHatSum = phase(complexSumHat);
thetaHatNorm = thetaHat - thetaHatSum;

% Calculate the shift theta
thetaSumTheory = atan2(sum(b),sum(a));
thetaHatSumTheory = atan2(sum(a),sum(b));
thetaShiftTheory = 90 - rad2deg(thetaSumTheory) - rad2deg(thetaHatSumTheory)

thetaShiftCalculate = rad2deg(thetaNorm(1) + thetaHatNorm(1))