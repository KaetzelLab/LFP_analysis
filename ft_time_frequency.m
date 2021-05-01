%% Input:
% indata: time points x channels x trials
% channels: 1 x 2 cell containing desired channels, e.g. 'PFCu', 'dHC'
%% Output:
% out: structure that contains structures with power, coherence, wPLI,
% Granger Causality (GC) and partial directed coherence (PDC)
% time-frequency spectra
%%
function [out] = ft_time_frequency(indata, channels)

data = {};
data.label = {'Ref'; 'PFCd'; 'PFCu'; 'dHC'; 'vHCu'; 'vHCd'; 'MD'}; % all channels
data.fsample = 1000; % sampling rate
Index1 = find(contains(data.label,channels(1))); % index of first channel
Index2 = find(contains(data.label,channels(2))); % index of second channel
idx = [];
for i = 1:size(indata,3) % loop through all trials and find trials where both channels are available (non-nan)
    if all(isnan(indata(:,Index1,i))) || all(isnan(indata(:,Index2,i)))
        idx = [idx,i];
    end
end
indata(:,:,idx) = []; % take only non-nan trials
for k = 1:size(indata,3) % import indata into FieldTrip structure
    data.time{k} = 1:7000;
    data.trial{k} = indata(:,:,k)';
end

cfg = [];
cfg.demean = 'yes';
data = ft_preprocessing(cfg, data);

cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = channels;
cfg.method       = 'wavelet';
cfg.foi          = 1:1:48;
cfg.toi          = 1:0.01:8;
cfg.width      = 3;
cfg.keeptrials = 'yes';
freq_data = ft_freqanalysis(cfg, data);

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = channels;
cfg.method       = 'wavelet';
cfg.foi          = 1:1:48;
cfg.toi          = 1:0.01:8;
cfg.width      = 3;
cfg.keeptrials = 'no';
pow_str = ft_freqanalysis(cfg, data);

%     figure
%     imagesc(10*log(squeeze(pow_str2.powspctrm(2,:,:))));
%     colorbar;colormap(jet)
%
%     figure; plot(pow_str.freq, 10*log(pow_str.powspctrm(1,:))); hold on; plot(pow_str.freq, 10*log(pow_str.powspctrm(2,:)))

cfg_conn = [];
cfg_conn.method = 'coh';
coh_str = ft_connectivityanalysis(cfg_conn, freq_data); 

%     figure
%     imagesc(squeeze(coh_str.cohspctrm(1,2,:,:)))
%
%     figure
%     plot(coh_str.freq, squeeze(coh_str.cohspctrm(1,2,:)))

cfg_conn = [];
cfg_conn.method = 'wpli_debiased';
wpli_str = ft_connectivityanalysis(cfg_conn, freq_data);

%     figure
%     imagesc(squeeze(wpli_str.wpli_debiasedspctrm(1,2,:,:)))

cfg_conn = [];
cfg_conn.method = 'granger';
granger_str = ft_connectivityanalysis(cfg_conn, freq_data); 

%     figure
%     imagesc(squeeze(granger_str.grangerspctrm(2,1,:,:)))

cfg_conn = [];
cfg_conn.method = 'pdc';
pdc_str = ft_connectivityanalysis(cfg_conn, freq_data); 

out.freq = pow_str.freq;
out.pow = pow_str.powspctrm;
out.coh = coh_str.cohspctrm;
out.wpli = wpli_str.wpli_debiasedspctrm;
out.GC = granger_str.grangerspctrm;
out.pdc = pdc_str.pdcspctrm;


figure(1)
subplot(4,2,1)
imagesc(squeeze(pow_str.powspctrm(2,:,:)))
title('Power CH2')
subplot(4,2,2)
imagesc(squeeze(pow_str.powspctrm(1,:,:)))
title('Power CH1u')
subplot(4,2,3)
imagesc(squeeze(wpli_str.wpli_debiasedspctrm(1,2,:,:)))
title('wPLI')
subplot(4,2,4)
imagesc(squeeze(coh_str.cohspctrm(1,2,:,:)))
title('Coherence')
subplot(4,2,5)
imagesc(squeeze(granger_str.grangerspctrm(1,2,:,:)))
title('CH1 -> CH2, GC')
subplot(4,2,6)
imagesc(squeeze(granger_str.grangerspctrm(2,1,:,:)))
title('CH2 -> CH1, GC')
subplot(4,2,7)
imagesc(squeeze(pdc_str.pdcspctrm(1,2,:,:)))
title('CH1 -> CH2, PDC')
subplot(4,2,8)
imagesc(squeeze(pdc_str.pdcspctrm(2,1,:,:)))
title('CH2 -> CH1, PDC')

