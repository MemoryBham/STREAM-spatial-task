% settings_retest

% nstim = 2*nstim;
% nsessions = 1;

distr_duration = 0;
retrievalTimeout = 0; %s


% number of times each stimulus is shown
nrep_fam = 0;
nrep_enc = 0; 
nrep_ret = 10;

% number of trials between two repetitions of the same stimulus
mindiff = 3;

% catch
nrep_Ctrials = 10; % number of catch trials for each stimulus;
mindiffC = 0;

% number of trials per block - there is a break after each block
ntrials_block_ret = 40;
if ntrials_block_ret < nstim
    warning('Number of trials per retrieval block is lower than the number of stimuli, the stimulus numbers will not be balanced!')
end
