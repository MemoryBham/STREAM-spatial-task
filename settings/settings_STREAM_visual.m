% settings_STREAM

% Initialize the random number generator
seed = sum(100*clock);
rng(seed);

%% paths

resultsFolder = './results/';

imageFolder = './stimuli/stimuli_dogs_birds_cars_planes/';

%% timing 

% reading time for questions during familiarization phase (seconds)
famTimeout = 2;

% How long showing each cue (in seconds)
cueTimeout = 1;

% How long showing each object (in seconds)
objectTimeout = 1.0;%0.5;

% Inter trial interval (in seconds)
postobjectTimeout = 1.0;%1.5;

% catch question duration
catchTimeout = 2;

% How long we show the feedback screen (in seconds)
feedbacktimeout = 1;

% How long (in seconds) is the minimum duration of the fixation point
fixationDuration = 0.5;

minBreakTime = 10; %s

retrievalTimeout = 0; %s

% distractor task length
distr_duration = 60; %s

%% trial numbers

nsessions = 1;

% number of times each stimulus is shown
nrep_fam = 25; % each question 5 times
nrep_enc = 0; 
nrep_ret = 0;

% number of trials between two repetitions of the same stimulus
mindiff = 3;

% drag-and-drop
nrep_DDtrials = 3; % number of drag-and-drop trials for each stimulus
mindiffDD = 1; % minimum number of trials, before and after a DD trials, that cannot be a DD trial

% catch
nrep_Ctrials = 10; % number of catch trials for each stimulus;
mindiffC = 0;

% number of trials per block - there is a break after each block
ntrials_block_enc = 16;
ntrials_block_ret = 16;
if ntrials_block_enc < nstim
    warning('Number of trials per encoding block is lower than the number of stimuli, the stimulus numbers will not be balanced!')
end
if ntrials_block_ret < nstim
    warning('Number of trials per retrieval block is lower than the number of stimuli, the stimulus numbers will not be balanced!')
end

%% catch questions

catch_types = {'exemplar', 'perc_1', 'perc_2', 'sem_1', 'sem_2'};

if mod(nrep_Ctrials, length(catch_types)) ~= 0
    warning('Catch frequency and number of question types is not balanced.\nSome catch types will appear more than others')
end
if mod(nrep_fam, length(catch_types)) ~= 0 
    warning('Familiarization frequency and number of question types is not balanced.\nSome questions will appear more than others')
end

catch_positions = '[(W/3) 1;(2*W/3) 3]'; % left, right
catch_y = 'H-(H/2-R)/2';%

%% keyboard inputs

% Keyboard setup
KbName('UnifyKeyNames');

% Response keys (optional; f or no subject response use empty list)
responseKeys = {'UpArrow','DownArrow','LeftArrow','RightArrow'};
stopKey = {'Q'};
instrKey = {'space'};
KbCheckList = [KbName(stopKey),KbName(instrKey)];%KbName('space'),KbName('ESCAPE'),

for i = 1:length(responseKeys)
    KbCheckList = [KbName(responseKeys{i}),KbCheckList];
end

RestrictKeysForKbCheck(KbCheckList);

%% some layout

fontsize = 40;

% for the presentation arena
R = 250; % radius of presentation circle, hardcoded on line 80
angle_offset = [0,pi/nstim];
dotsize = 5; % size for dots indicating locations

% text offset
instroffset = 500; % pixels
labeloffset = -25;

% image size in pixels
imageSize_o = [300,300];

% for breaks between blocks
barWidth = 500;
barHeight = 50;

%% colors

% Background color: choose a number from 0 (black) to 255 (white)
backgroundColor = 200;

% Text color: choose a number from 0 (black) to 255 (white)
textColor = 0;

% fixation cross
barColor = [0 0 0]; %Controlling fix cross color when there is no feedback

% cue indicator color
barColorCue = [0,0,0];%[0 50 155];

white = 255*[1,1,1]; %WhiteIndex(window1);
black = [0,0,0];%BlackIndex(window1);

% colors for feedback during drag-and-drop encoding trials
barColorDD = [0.55,0.55,0.55]*255;

% colors for correct and incorrect answers
barColor_Correct = [0.1 0.5 0.2]*255;%[0 1 0]*255; % green
barColor_Incorrect = [190 30 45];%[1 0 0]*255; % red
