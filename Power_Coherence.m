%% Calculates Power, Coherence and phase angle using multitapers
% Input: indata (Nsamples x Nchannels)
% Output: Power and Coherence spectrum (see mtspectrumc and coherencyc for
% further information)
%
% Reference: Chronux toolbox: Bokil et al., J Neurosci Methods 2010
% (http://chronux.org/)

% Daniel Strahnen, 2020

%% params for spectral analysis
params.Fs = 1000; % sampling rate
params.fpass = [0 48]; %frequency range for which to analyse
params.pad = -1; % no padding
params.trialave = 1; % 1 = average plot 0 = plots individual traces
params.tapers = [0.2,600,20]; % define tapers
params.err = [2 0.05]; % Jackknife error bars (required for the coherence plots)

[S,f,Serr] = mtspectrumc(indata(:,1),params); % Power spectrum for Channel 1
S = 10*log(S);

[C,phi,S12,S1,S2,f,confC,phistd,Cerr] = coherencyc(indata(:,3),indata(:,4),params); % Coherence and phase angle phi between Channel 3 and Channel 5
