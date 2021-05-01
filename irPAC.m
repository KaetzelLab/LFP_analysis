%% Calculates phase-amplitude coupling using the Modulation Index
% Input: data1 - data array from Region 1
%        data2 - data array from Region2
% Output: out - structure that contains the MI (for further information see
%         get_mi)
% References: Tort et al., PNAS 2008, 2009
%             Voloh et al., PNAS 2015
%             (MatLab code downloaded from: https://github.com/bvoloh/cfc_analysis)
%             eegfilt: EEGLAB: Delorme & Makeig, J Neurosci Methods 2004, https://sccn.ucsd.edu/~arno/eeglab/auto/eegfilt.html

srate = 1000;
lo_freq_phase = 5;
hi_freq_phase = 12;
lo_freq_amp = 30;
hi_freq_amp = 48;
surrogate_runs = 100;

phasedata = eegfilt(data1', srate, lo_freq_phase, []);
phasedata = eegfilt(phasedata, srate, [], hi_freq_phase);
phasedata = angle(hilbert(phasedata)); % phase

ampdata = eegfilt(data2', srate, lo_freq_amp, []);
ampdata = eegfilt(ampdata, srate, [], hi_freq_amp);
ampdata = abs(hilbert(ampdata)); % amplitude

%% MI
%inputs
nbins=18;
randtype=2; % slices

out = get_mi(phasedata',ampdata',nbins,surrogate_runs,randtype);