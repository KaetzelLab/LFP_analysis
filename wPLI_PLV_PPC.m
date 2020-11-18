%% Calculates weighted phase lag index (wPLI), phase-locking value (PLV) and 
%% pairwise phase consistency (PPC) based on the configurations below
% Input: indata (Nsamples x Nchannels)
% Output: wpli_str, plv_str and ppc_str contain respective spectra
%
% References: 
% FieldTrip toolbox: Oostenveld et al., Comput Intell Neurosci 2010 (https://www.fieldtriptoolbox.org/)
% wPLI: Vinck et al., Neuroimage 2011
% PLV: Lachaux et al., Hum Brain Mapp 1999
% PPC: Vinck et al., Neuroimage 2010

% Daniel Strahnen, 2020

%% Convert indata to data structure according to FieldTrip toolbox
data = {};
data.label = {'Ref'; 'PFCd'; 'PFCu'; 'dHC'; 'vHCu'; 'vHCd'; 'MD'}; % Recorded Channels
data.fsample = 1000;
j = 1;
for k = 1:1
    data.time{k} = j:(k*600000);
    data.trial{k} = indata((j:(k*600000)),:)';
    j = k*600000 + 1;
end
%% Change data length to 1s with 50% overlap
cfg.length    =  1;
cfg.overlap   =  0.5;
[data] = ft_redefinetrial(cfg, data);
%% Freqanalysis to obtain the cross spectral density
cfg_freq = [];
cfg_freq.method = 'mtmfft';
cfg_freq.output = 'fourier';
cfg_freq.keeptrials = 'yes';
cfg_freq.taper = 'hann';
cfg_freq.tapsmofrq = 0.5;
cfg_freq.foilim = [0 48];
cfg_freq.pad = 'nextpow2';
freq_data = ft_freqanalysis(cfg_freq, data);
%% debiased wPLI analysis
cfg_conn = [];
cfg_conn.method = 'wpli_debiased';
wpli_str = ft_connectivityanalysis(cfg_conn, freq_data);
wpli_str = ft_checkdata(wpli_str, 'cmbrepresentation', 'full','datatype','freq');
%% PLV
cfg_conn = [];
cfg_conn.method = 'plv';
plv_str = ft_connectivityanalysis(cfg_conn, freq_data);
plv_str = ft_checkdata(plv_str, 'cmbrepresentation', 'full','datatype','freq');
%% PPC
cfg_conn = [];
cfg_conn.method = 'ppc';
ppc_str = ft_connectivityanalysis(cfg_conn, freq_data);
ppc_str = ft_checkdata(ppc_str, 'cmbrepresentation', 'full','datatype','freq');