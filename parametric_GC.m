%% Calculates parametric Granger Causality using the MVGC-toolbox
% Input: data - Nchannels x Nsamples x Ntrials, data sampled down to 250Hz
% Output: Matrix containing GC spectra and integrated GC values over
% certain frequency bands
% Reference: Barnett et al., J Neurosci Methods 2014 (https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html)
% Code is largely based on the provided demo version

% Daniel Strahnen, 2020

morder = 27; % set model order

ntrials   = size(data,3);     % number of trials
nobs      = size(data,2) ;   % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)

fs        = 250;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

rng_seed(seed);

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)
ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(data,morder,regmode);
ptoc;

assert(~isbad(A),'VAR estimation failed');

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>
ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;
var_info(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

assert(~isbad(F,false),'GC calculation failed');
%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

ptic('\n*** autocov_to_spwcgc... ');
f = autocov_to_spwcgc(G,fres);
ptoc;

assert(~isbad(f,false),'spectral GC calculation failed');

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
mad = maxabs(F-Fint);
madthreshold = 1e-5;
if mad < madthreshold
fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
else
    fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
end

%% Aggregated ("band-limited") GC in specific key frequency bands
Nyquist = fs / 2;
freq_range = [1/Nyquist, 4/Nyquist]; % freq_range is a vector that can reach from 0 to 1 (Nyqvist)
GC_delta = smvgc_to_mvgc(f, freq_range);
freq_range = [5/Nyquist, 12/Nyquist]; % freq_range is a vector that can reach from 0 to 1 (Nyqvist)
GC_theta = smvgc_to_mvgc(f, freq_range);
freq_range = [15/Nyquist, 30/Nyquist]; % freq_range is a vector that can reach from 0 to 1 (Nyqvist)
GC_beta = smvgc_to_mvgc(f, freq_range);
freq_range = [30/Nyquist, 48/Nyquist]; % freq_range is a vector that can reach from 0 to 1 (Nyqvist)
GC_gamma = smvgc_to_mvgc(f, freq_range);

%% Permutation test
alpha     = 0.05;   % significance level for all statistical tests
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
bsize     = [];     % permutation test block size: empty for automatic (uses model order)
nperms    = 100;    % number of permutations for permutation test
ptic('\n*** tsdata_to_mvgc_pwc_permtest\n');
FNULL = permtest_tsdata_to_pwcgc(data,morder,bsize,nperms,regmode,acmaxlags);
ptoc('*** tsdata_to_mvgc_pwc_permtest took ',[],1);

pval_p = empirical_pval(GC_delta,FNULL);
sig(:,:,1)  = significance(pval_p,alpha,mhtc);
pval_p = empirical_pval(GC_theta,FNULL);
sig(:,:,2)  = significance(pval_p,alpha,mhtc);
pval_p = empirical_pval(GC_beta,FNULL);
sig(:,:,3)  = significance(pval_p,alpha,mhtc);
pval_p = empirical_pval(GC_gamma,FNULL);
sig(:,:,4)  = significance(pval_p,alpha,mhtc);

%% Save the data into a matrix
f = permute(f, [3,1,2]);
freq = linspace(0,Nyquist,size(f,1)); % freq goes until Nyqvist, so 125Hz for sample rate of 250Hz
