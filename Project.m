


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
% 8–30 Hz to remove the molestation information. The design specifications for 
% Butterworth band-pass filter are: 
% the stop-band cutoff frequency is 5 Hz and 33 Hz respectively, 
% stop-band attenuation is 50 dB?
% and pass-band attenuation is 0.5 dB?
sample_rate = EEG.srate;
pass_band = [8 30];         % they dont use that in the paper?
stop_band = [5 33];         % stop-band cutoff
order = 6;

% butterworth band-pass filter
[b, a] = butter(order, stop_band / (sample_rate / 2), 'stop');
% ref: https://www.mathworks.com/help/signal/ref/butter.html

% apply the filter to data
EEG.data = filter(b, a, EEG.data, [], 2);


%% 3. Feature extraction - try implement FAWT as in paper
% In the current paper, parameters of FAWT are p = 5, q = 6, r = 1,
% s = 2 and beta = (0.8r)/s [29–31], and the level of decomposition is
% kept at J = 5
p = 5; 
q = 6; 
r = 1; 
s = 2; 
beta = (0.8 * r) / s;
levels_of_dec = 5;

function sub_bands = fawt(signal, levels_of_dec, p, q, r, s, beta)
    % paper pg.3 - ch.2.3
    % Using iterative filter bank, the raw signal is decomposed into 5th level

    % low pass filter
        % H(w)
    % high pass filter
        % G(w)
end

% FAWT on each channel
num_channels = size(EEG.data, 1);
num_samples = size(EEG.data, 2);
% ref: https://www.mathworks.com/help/matlab/ref/double.size.html

% F is a 60-dimensions feature vector, but we construct matrix, because of
% the channels
features = zeros(num_channels, 60); 

for ch = 1:num_channels
    signal = EEG.data(ch, :);
    sub_bands = fawt(signal, levels_of_dec, p, q, r, s, beta);
    % !important!
    % the sub-band signals reconstructed by FAWT are denoted as SB
   
    % get features from each sub-band
    feature_vector = [];
    for sb = 1:levels_of_dec
        % get sb signal
        sb_signal = sub_bands(sb, :); 
        % calculate mean energy
        T_E = mean(sb_signal.^2);
        % calc      abs mean value
        T_aav = mean(abs(sb_signal));
        % calc      standard deviation
        T_std = std(sb_signal);
        % calc      voltality idx
        T_vi = sum(abs(diff(sb_signal))) / (num_samples - 1);
        % calc      center freq
        P_cf = calculate_center_frequency(sb_signal, EEG.srate);
        
        % append all five parameters features
        T_cov = 10;
        % what is T_cov?
        feature_vector = [feature_vector, T_E, T_std, T_cov, P_cf];
    end
    features(ch, :) = feature_vector;
end

% construct a feature vector - look at pg.2 (5) in paper
% F = [TE Tstd Tcov Pcf ]
% should get F 60 dimension vector

%% 4. feature reduction of vector F - MDS (razišči kako deluje)
% Construction of distance matrix D. For each pair of
% vectors xi and xj (1 ≤ j ≤ N), the distance matrix D = (di,j)
% N×N is measured using Euclidean distance di,j
% - create matrix D and use provided mdsclae function
D = pdist(features, 'euclidean');
% ref: https://www.mathworks.com/help/stats/pdist.html

% For each subject a traversal searching for d from 2 to 60 is performed on MDS
p = 5; % change when testing
[reduced_features, ~] = mdscale(D, p);
% ref: https://www.mathworks.com/help/stats/mdscale.html

%% 5. Classification with LDA and others classificators 
% what we get out of mdscale function we use as an input to LDA
% provided in matlab: MdlLinear = fitcdiscr(meas,species);

%% 6.Compare classificators between each other
% 