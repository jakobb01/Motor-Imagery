


%% 1. LOAD DATA
eeglab
filepath= 'C:\Users\jakob\Documents\MATLAB\eeg-motor-movementimagery-dataset-1.0.0\files\S001\';
filename='S001R01.edf'; %filename={'S001R01.edf', 'S001R02.edf'};

EEG = pop_biosig([filepath filesep filename]);
EEG = pop_chanedit(EEG, 'load',{['C:/Users/jakob/Documents/MATLAB/eeg-motor-movementimagery-dataset-1.0.0/files/BCI2000.locs'] 'filetype' 'autodetect'});
EEG = pop_reref( EEG, []);

% load events, as EEG.event is empty


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

    % use formulas provided in paper to create fawt.
    % low pass filter
        % H(w)
    % high pass filter
        % G(w)
    current_signal = signal;
    sub_bands = cell(levels_of_dec + 1, 1); % init sub_bands
    
    for level = 1:levels_of_dec
        % get filters
        [H, G] = fawt_filters(p, q, r, s, beta, level);
        % apply filters
        low_pass = filter(H, 1, current_signal);
        high_pass = filter(G, 1, current_signal);
        % according to schematic on the ref paper
        sub_bands{level} = high_pass;
        current_signal = low_pass;
    end
    % store the last low pass one as a sub_band
    sub_bands{levels_of_dec + 1} = current_signal;
        
end

function [H, G] = fawt_filters(p, q, r, s, beta, level)
    % filter parameters based of the paper
    wp = ((1 - beta) * pi + eps) / p;
    ws = pi / q;
    w0 = ((1 - beta) * pi + eps) / r;
    w1 = (p * pi) / (q * r);
    w2 = (pi - eps) / r;
    w3 = (pi + eps) / r;
    function O = help_fun(w)
        O = sqrt(2 - cos(w)) .* (1 + cos(w)) / 2;
    end

    w = linspace(0, 2 * pi, 1024);

    % low pass filter
    H = zeros(1, length(w));
    for i = 1:length(w)
        if (abs(w(i)) < wp)
            H(i) = sqrt(p*q);
        elseif (wp <= w(i) && w(i) <= ws)
            H(i) = sqrt(p*q) * help_fun(((w(i) - wp) / (ws - wp)));
        elseif (-ws <= w(i) && w(i) <= -wp)
            H(i) = sqrt(p*q) * help_fun(((pi - w(i) + wp) / (ws - wp)));
        elseif (abs(w(i)) > ws)
            H(i) = 0;
        else 
            H(i) = 0;
            disp('H found no condition');
        end
    end
    % high pass filter
    G = zeros(1, length(w));
    for i = 1:length(w)
        if (w0 <= w(i) && w(i) <= w1)
            G(i) = sqrt(2*r*s) * help_fun(((pi - w(i) - w0) / (w1 - w0)));
        elseif (w1 < w(i) && w(i) < w2)
            G(i) = sqrt(2*r*s);
        elseif (w2 <= w(i) && w(i) <= w3)
            G(i) = sqrt(2*r*s) * help_fun(((w(i) - w2) / (w3 - w2)));
        elseif ((0 <= w(i) && w(i) < w0) || (w3 < w(i) && w(i) <= 2 * pi))
            G(i) = 0;
        else 
            G(i) = 0;
            disp('G found no condition');
        end
    end
    % freq response from transfer function with inverse FFT
    H = ifft(ifftshift(H));
    G = ifft(ifftshift(G));
    % old ref: https://www.mathworks.com/help/signal/ref/freqz.html
end

function P_cf = center_frequency(sb_signal, sample_rate)
    % as they describe in the paper
    N = length(sb_signal);
    freq = (0:N-1) * (sample_rate / N); 
    pow = abs(fft(sb_signal)).^2; % Fourier transform seperately
    P_cf = sum(freq .* pow) / sum(pow);
end

% FAWT on each channel
num_channels = size(EEG.data, 1);
num_samples = size(EEG.data, 2);
% ref: https://www.mathworks.com/help/matlab/ref/double.size.html

% F is a 60-dimensions feature vector, but we construct matrix, because of
% the channels
features = zeros(num_channels, (levels_of_dec + 1) * 4);

for ch = 1:num_channels
    signal = EEG.data(ch, :);
    sub_bands_cell = fawt(signal, levels_of_dec, p, q, r, s, beta);
    % !important!
    % the sub-band signals reconstructed by FAWT are denoted as SB
   
    % get features from each sub-band
    feature_vector = [];
    for sb_idx = 1:length(sub_bands_cell )
        % get sb signal
        sb_signal = sub_bands_cell {sb_idx};
        if ~isempty(sb_signal) % Ensure sub-band signal is not empty
            % calculate mean energy
            T_E = mean(sb_signal.^2);
            % calc      abs mean value
            % T_aav = mean(abs(sb_signal));
            % calc      standard deviation
            T_std = std(sb_signal);
            % calc      voltality idx
            T_vi = sum(abs(diff(sb_signal))) / (num_samples - 1);
            % calc      center freq
            P_cf = center_frequency(sb_signal, sample_rate);
            
            % append all five parameters features
            % T_cov = T_vi. Mention in the paper -> pg.2 (introduction)
            feature_vector = [feature_vector, T_E, T_std, T_vi, P_cf];
        end
    end
    features(ch, :) = feature_vector;
end

% construct a feature vector - look at pg.2 (5) in paper
% F = [TE Tstd Tcov Pcf ]
% should get F 60 dimension vector

%% 4. feature reduction of vector F - MDS
% Construction of distance matrix D. For each pair of
% vectors xi and xj (1 ≤ j ≤ N), the distance matrix D = (di,j)
% N×N is measured using Euclidean distance di,j
% - create matrix Ds
D = pdist(real(features), 'euclidean');
% ref: https://www.mathworks.com/help/stats/pdist.html

% For each subject a traversal searching for d  from 2 to 60 is performed on MDS
p = 5; % change when testing
[reduced_features, ~] = mdscale(D, p);
% ref: https://www.mathworks.com/help/stats/mdscale.html

%% 5. Classification with LDA and others classificators 
% what we get out of mdscale function we use as an input to LDA
% provided in matlab: MdlLinear = fitcdiscr(meas,species);

% create labels


% split the data into training and testing sets


train_indices = 0;
test_indices = 0;

train_features = reduced_features(train_indices, :);
test_features = reduced_features(test_indices, :);
train_labels = labels(train_indices);
test_labels = labels(test_indices);

% train LDA classifier
lda_model = fitcdiscr(train_features, train_labels);

% calculate accuracy
predicted_labels = predict(lda_model, test_features);
accuracy = sum(predicted_labels == test_labels) / length(test_labels) * 100;
fprintf('Classification Accuracy: %.2f%%\n', accuracy);

% Visualize


%% 6.Compare classificators between each other
% 