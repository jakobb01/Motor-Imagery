


%% 1. LOAD DATA
eeglab
filepath= 'C:\Users\jakob\Documents\MATLAB\eeg-motor-movementimagery-dataset-1.0.0\files\S001\';
filename='S001R01.edf'; %filename={'S001R01.edf', 'S001R02.edf'};

EEG = pop_biosig([filepath filesep filename]);
EEG = pop_chanedit(EEG, 'load',{['C:/Users/jakob/Documents/MATLAB/eeg-motor-movementimagery-dataset-1.0.0/files/BCI2000.locs'] 'filetype' 'autodetect'});
EEG = pop_reref( EEG, []);


%% 2. PREPROCESSING
% paper description:
% mainly occurred
% in Alpha (8–13 Hz) and Beta (14–28 Hz) rhythms [24,25], a six order
% Butterworth band-pass filter is applied on the raw EEG signals in
% 8–30 Hz to remove the molestation information. The design specifications for Butterworth band-pass filter are: the stop-band cutoff
% frequency is 5 Hz and 33 Hz respectively, stop-band attenuation is
% 50 dB, and pass-band attenuation is 0.5 dB


%% 3. Feature extraction - try implement FAWT as in paper
% In the current paper, parameters of FAWT are p = 5, q = 6, r = 1,
% s = 2 and ˇ = (0.8r)/s [29–31], and the level of decomposition is
% kept at J = 5


% !important!
% the sub-band signals reconstructed by FAWT are denoted as SB

% construct a feature vector - look at pg.2 (5) in paper
% F = [TE Tstd Tcov Pcf ]
% should get F 60 dimension vector

%% 4. feature reduction of vector F - MDS (razišči kako deluje)
% Construction of distance matrix D. For each pair of
% vectors xi and xj (1 ≤ j ≤ N), the distance matrix D = (di,j)
% N×N is measured using Euclidean distance di,j
% - create matrix D and use provided mdsclae function


%% 5. Classification with LDA or others classificators
% what we get out of mdscale function we use as an input to LDA
% provided in matlab: MdlLinear = fitcdiscr(meas,species);