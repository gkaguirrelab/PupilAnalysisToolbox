%% MacroRunAllAnalyses.m
% Runs analysis of all pupillometry protocols
%
% 10/12/13  ms      Commented.

%% Run analysis of TTF4D protocol
TTF4D_100secWindow_Analysis_Params;         
TTF4D_100secWindow_HarmonicAnalysis_Params; % Harmonic analysis

%% Run analysis of OPN4 protocols
OPN4_40secWindow_05Hz_Analysis_Params;      % 0.5 Hz
OPN4_40secWindow_005Hz_Analysis_Params;     % 0.05 Hz

%% Run analysis of LumVar protocols
LumVarND10_100secWindow_Analysis_Params;            % ND1.0
LumVarND10_100secWindow_HarmonicAnalysis_Params;    % + Harmonic
LumVarND20_100secWindow_Analysis_Params;            % ND2.0
LumVarND20_100secWindow_HarmonicAnalysis_Params;    % + Harmonic

%% Run analysis of MFC protocol
MFC_100secWindow_Analysis_Params;

%% Run analysis of RISC protocols (robust S)
RISC_100secWindow_Analysis_Params;
RISC_100secWindow_HarmonicAnalysis_Params;

%% Run analysis of SMelPAT protocol
SMelPAT_100secWindow_Analysis_Params;

%% Run analysis of SPP protocol
SPP_25secWindow_Analysis_Params;

%% Run analysis of envelope pupil stuff
Envelope_40secWindow_05Hz_Analysis_Params;

