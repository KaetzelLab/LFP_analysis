%% Calculates non-parametric GC, directed transfer function (DTF) and partial directed coherence (PDC)
%% Input: data - Nsamples x Nchannels x Ntrials, sampled down to 250Hz
%% Output: structures containing respective spectral information
%% References:
%  Non-parametric Granger Causality: Dhamala et al., Phys Rev Lett 2008 and
%  Dhamala et al., Neuroimage 2008
%  DTF: Kaminski & Blinowska, Biol Cybern 1991
%  PDC: Baccala et al., Biol Cybern 2001
%  MatLab functions: FieldTrip toolbox: Oostenveld et al., Comput Intell Neurosci 2010 (https://www.fieldtriptoolbox.org/)

%  Daniel Strahnen, 2020

%% Convert indata to data structure according to FieldTrip toolbox
data_raw = {};
data_raw.label = {'CH1'; 'CH2'}; % 2 example channels
data_raw.fsample = 250;
for k = 1:size(data,3)
    data_raw.time{k} = 1:size(data,1);
    data_raw.trial{k} = data(:,:,k)';
end
%% Change data length to 1s with 50% overlap
cfg.length    =  1;
cfg.overlap   =  0.5;
[data_raw] = ft_redefinetrial(cfg, data_raw);
%% Freqanalysis to obtain the cross spectral density
cfg_freq = [];
cfg_freq.method = 'mtmfft';
cfg_freq.output = 'fourier';
cfg_freq.keeptrials = 'yes';
cfg_freq.taper = 'dpss';
cfg_freq.tapsmofrq = 2;
cfg_freq.foilim = [0 48];
cfg_freq.pad = 'nextpow2';
freq_data = ft_freqanalysis(cfg_freq, data_raw);

%% directionality analysis
cfg_conn = [];
cfg_conn.method = 'granger';
granger_str = ft_connectivityanalysis(cfg_conn, freq_data);

cfg_conn = [];
cfg_conn.method = 'pdc';
pdc_str = ft_connectivityanalysis(cfg_conn, freq_data);

cfg_conn = [];
cfg_conn.method = 'dtf';
dtf_str = ft_connectivityanalysis(cfg_conn, freq_data);