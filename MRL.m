%% Calculates MUA-Phase Locking
%% Inputs:
% indata1: data array for MUA, e.g. PFC
% indata2: data array for LFP, e.g. vHipp
%% Output:
% results: Array that contains the mean angle, MRL, p-value (Rayleigh's
% test), z-value and lead or lag with maximum MRL
%% References:
% Methodology: Siapas et al., Neuron 2005; Sigurdsson et al., Nature 2010
% MatLab functions: CircStat toolbox: Berens, J Stat Softw 2009
% (https://github.com/circstat/circstat-matlab)
% eegfilt: EEGLAB: Delorme % Makeig, J Neurosci Methods 2004
% (https://sccn.ucsd.edu/~arno/eeglab/auto/eegfilt.html)

%  Daniel Strahnen, 2020

%% Define frequency ranger and filter data
hi_cutoff = 800;
low_freq = 5;
high_freq = 12;
srate = 10000;
hp = highpass(indata1,hi_cutoff,srate); 
theta = eegfilt(Output_data_raw(:,2)',srate,low_freq,high_freq);
%% Interpolate between peaks because of waveform asymmetry (Buzsaki et al., Electroencephalogr Clin Neurophysiol 1985; McClain et al., PNAS 2019)
locmin = find(islocalmin(theta)); % find troughs of theta
thetaphase = zeros(1,length(theta));
for i = 1:length(locmin) % interpolate phase between troughs
    if (i+1) <= length(locmin)
        thetaphase(locmin(i):locmin(i+1)-1) = linspace(-pi,pi,locmin(i+1)-locmin(i));
    end
end
for i = 1:length(thetaphase) % shift thetaphase so that a phase of zero refers to the trough of the theta cycle (Adhikari et al., Neuron 2010)
    if thetaphase(i) > 0
        thetaphase_shifted(i) = thetaphase(i) - pi;
    else
        thetaphase_shifted(i) = thetaphase(i) + pi;
    end
end
%% Determine threshold for theta amplitude
theta_amp1 = abs(hilbert(theta));
theta_amp1 = theta_amp1';
theta_amp = theta_amp1-mean(theta_amp1);
theta_threshold = 0.25*std(theta_amp);
%% Spikes with changable std
threshold = 3.5*std(hp,0,1); % threshold for spike detection
for i = 1:size(hp,1)
    if hp(i) > threshold
        B(i,1) = 1;
    else
        B(i,1) = 0;
    end
end
downsample_factor = srate / 1000; % to 1kHz
%% Remove spikes that occur within 1 ms (10 data points) of each other or are longer than 2 ms (20 data points)
k = 1;
k_start_previous = 1;
while k <= size(B,1)
    while (k <= size(B,1)) && (B(k) == 0)
        k = k + 1;
    end
    k_start = k;
    k_end = k_start;
    while (k_end <= size(B,1)) && (B(k_end) == 1)
        k_end = k_end + 1;
    end
    if (k_end - k_start) > 20 % check if spike is longer than 2ms
        B(k_start:k_end) = 0;
    else
        if (k_start - k_start_previous) < 10 % check if spikes occured within 1ms
            B(k_start:k_end) = 0;
            B(k_start_previous) = 0;
        else
            B(k_start+1:k_end) = 0; % keep only first timepoint of the spike
        end
    end
    k_start_previous = k_start;
    k = k_end + 1;
end
%% Take only spikes when theta was above threshold
for i = 1:length(B)
    if B(i) == 1 && (i-0.1*srate) > 1 && (i+0.1*srate) < length(B)
        j = -0.1*srate;
        while j <= 0.1*srate
            if theta_amp(i+j) < theta_threshold
                B(i) = 0;
                j = 0.1*srate + 1;
            else
                j = j + 1;
            end
        end
    end
end
B_idx = find(B); % indices of spikes
%% Determine phase of each spike
phases = [];
for i = 1:length(B_idx)
    phases(i) = thetaphase_shifted(B_idx(i))+pi;
end
phases_rand = phases(randperm(length(phases))); % take random phases
numangles = 1000; % number of spikes to use

x = length(phases);

if length(phases) < numangles % return NaNs if criteria are not met
    results = nan(1,5);
    phases = nan(1,numangles);
    R_offset = nan(1,51);
else
    phases = phases(1,1:numangles);
    
    alpha_bar = circ_mean(phases'); % mean angle
    R = circ_r(phases'); % MRL
    [pval, z] = circ_rtest(phases'); % Rayleigh's test
    
    results = [R, alpha_bar, pval, z];
    %% Directionality by shifting the spikes relative to the thetaphase and get maximum
    p_offset_all = [];
    R_offset = [];
    lags = -0.1*srate:0.004*srate:0.1*srate;
    for i = -0.1*srate:0.004*srate:0.1*srate
        y = [];
        for j = 1:size(B,1)
            if (j+i < size(B,1)) && j+i > 0 && B(j) == 1
                y = [y, thetaphase_shifted(j+i)];
            end
        end
        [p_offset, ~] = circ_rtest(y');
        p_offset_all = [p_offset_all, p_offset];
        R_offset = [R_offset, circ_r(y')];
    end
    if ~isnan(sum(R_offset)) && ~isnan(sum(p_offset_all))
        R_offset = (R_offset - min(R_offset)) / (max(R_offset) - min(R_offset)); % normalize between 0 and 1
        g=find(R_offset==max(R_offset)); % identifies index where the MRL peaks
        R_offsetpeak=lags(g); % identifies the lag at which the MRL peaks
        results = [results, R_offsetpeak / downsample_factor]; % convert to ms
    else
        results = [results, nan];
    end
end
