function run_STREAMspatial()
clear all

DPIswitch = 0; % 0 for normal screens

%% the preliminaries

addpath('./settings/')
addpath('./functions/')

ExpStart = GetSecs;
today = date;
now = datestr(datetime('now'), 'dd-mmm-yyyy_HH.MM');

%% get settings

% Open a Prompt Window
prompt = {'subject ID', 'Number of stimuli (4,8,16)',...
    'Stimulus selection (auto/rand/[numb])',...
    'Session # (1, 2, 3 or 4)','Show instructions (yes/no)?',...
    'Trigger type (ttl, serial, labjack or none)', 'Run as practice round (yes/no)?',...
    'Task type (behavioral/standard/visual)'};

defaults = {'subXX','8','rand','1','yes','ttl','no','standard'};

answer = inputdlg(prompt, 'Experimental Setup Information', 1, defaults);
[subID, nstim, switch_selection, sessionID, switch_instructions, trg_type, practice, task_type] = deal(answer{:});

nstim = str2double(nstim);
sessionID = str2double(sessionID);

% get others settings
switch_behavior = false;
switch_standard = false;
switch_visual = false;
if sum(strcmp(practice, {'yes','y','Y','YES'})) % practice round
    switch_selection = 'rand';
    subID = strcat(subID, '_practice');
    if sum(strcmp(task_type, {'behavior','behaviour','behavioral','beh','BEH','Behavior','Behaviour'}))
        task_type = 'behavior';
        subID = strcat(subID, '_behavior');
        settings_STREAM_behavior_practice
        switch_behavior = true;
    elseif sum(strcmp(task_type, {'visual','vis','VIS','Visual'}))
        task_type = 'visual';
        subID = strcat(subID, '_visual');
        settings_STREAM_visual_practice
        switch_visual = true;
    else
        task_type = 'standard';
        settings_STREAM_practice
        switch_standard = true;
    end
    if nstim>4
        warning('Number of stimuli is higher than allowed for a practice round. Reducing nstim to 4.')
        nstim = 4;
    end
elseif sum(strcmp(practice, {'no','n','N','NO'})) % real session
    if sum(strcmp(task_type, {'behavior','behaviour','behavioral','beh','BEH','Behavior','Behaviour'}))
        if sessionID ~= 3 && sessionID ~= 4
            subID = strcat(subID, '_behavior');
        end
        task_type = 'behavior';
        settings_STREAM_behavior
        switch_behavior = true;
    elseif sum(strcmp(task_type, {'visual','vis','VIS','Visual'}))
        task_type = 'visual';
        if sessionID ~= 3 && sessionID ~= 4
            subID = strcat(subID, '_visual');
        end
        settings_STREAM_visual
        switch_visual = true;
    else
        task_type = 'standard';
        settings_STREAM
        switch_standard = true;
    end
    if sessionID == 3 || sessionID == 4
        settings_retest
    end
end

save(['./settings/settings_',subID,'_',now,'.mat'])

%% look for existing participant info

if exist(['.\results\resultfile_',subID,'.txt'],'file')
    % load existing log
    old_log = tdfread([resultsFolder '/resultfile_' subID '.txt']);
    
    if length(old_log.sessionID)==0
        fclose('all');
        error('Empty results file found, delete this file and start again.')
    end
    
    % check session ID
    if ischar(old_log.sessionID(end,1))
        old_session = str2double(old_log.sessionID(end,1));
    else
        old_session = old_log.sessionID(end,1);
    end
    if isempty(old_log.sessionID) || strcmp(old_log.sessionID(end,1),'s') || old_session~=sessionID
        cont = input('Data found for this participant, but not for this session. Continue with old stimulus-set? (yes/no/quit)\n','s');
        if sum(strcmp(cont,{'y','yes','Y', 'YES'}))
            new_participant = 0;
            tr_fam = 1;
            tr_enc = 1;
            tr_ret = 1;
        elseif sum(strcmp(cont,{'n','no','N', 'NO'}))
            new_participant = 1;
        elseif sum(strcmp(cont,{'q','stop','quit'}))
            close all
            sca
            return
        else
            error('Unknown input')
        end
        
    elseif ~strcmp(old_log.sessionID(end,1),'s') && old_session==sessionID
        %     if old_log.sessionID(end,1)==sessionID
        cont = input('Existing session found for this participant. Continue with the old session? (yes/no/quit)\n','s');
        if sum(strcmp(cont,{'y','yes','Y', 'YES'}))
            new_participant = 0;
            % find last trial for encoding and retrieval
            if ischar(old_log.trial_id(end))
                old_trial = str2double(old_log.trial_id(end,1));
            else
                old_trial = old_log.trial_id(end);
            end
            if strcmp('familiar',old_log.block_state(end,1:8))
                tr_fam = old_trial+1;
                tr_enc = 1;
                tr_ret = 1;
            elseif strcmp('encoding',old_log.block_state(end,1:8))
                tr_fam = nrep_fam*nstim+1;
                tr_enc = old_trial+1;%str2double(old_log.trial_id(end,1:5))+1;
                tr_ret = 1;
            elseif strcmp('retrieva',old_log.block_state(end,1:8))
                tr_fam = nrep_fam*nstim+1;
                tr_enc = nrep_enc*nstim+1;
                tr_ret = old_trial+1;
            end
        elseif sum(strcmp(cont,{'n','no','N', 'NO'}))
            new_participant = 1;
        elseif sum(strcmp(cont,{'q','stop','quit'}))
            close all
            sca
            return
        else
            error('Unknown input')
        end
    else
        new_participant = 1;
    end
else
    new_participant = 1;
end

if new_participant == 1
    tr_fam = 1;
    tr_enc = 1;
    tr_ret = 1;
end

%% Screen setup

clear screen
Screen('Preference', 'SkipSyncTests', 1);  %%% IMPORTANT!! REMOVE THIS LINE
whichScreen = max(Screen('Screens'));
[window1, rect] = Screen('Openwindow',whichScreen,backgroundColor,[],[],2);
slack = Screen('GetFlipInterval', window1)/2;
hz=Screen('NominalFrameRate', window1);

[xCenter, yCenter] = RectCenter(rect);
W=rect(RectRight); % screen width
H=rect(RectBottom); % screen

catch_positions = eval(catch_positions);
catch_y = eval(catch_y);

Screen(window1,'FillRect',backgroundColor);
Screen('Flip', window1);

%% get stimulus info

if new_participant
    % get image info from
    stimuli_info = readtable(strcat(imageFolder, '/stimuli_info.txt'), 'Delimiter', 'tab');
    
    if sum(strcmp(switch_selection, {'rand', 'random'}))
        % assign stimuli to two sessions
        stimuli_info = create_stimulus_split(stimuli_info,nsessions, nstim);
        
        % find locations for both sessions
        stimuli_info.positions = zeros(size(stimuli_info,1),2);
        for s = 1:nsessions
            [cue_positions,cue_positions_rand] = create_circular_array_angleoffset...
                (W/2,H/2,R,angle_offset(s),stimuli_info(stimuli_info.session == s,:));
            stimuli_info.positions(stimuli_info.session==s,:) = cue_positions_rand;
        end
        selectionID = 0;
    else
        % follow predefined selection
        load([imageFolder, '/stimuli_locations.mat'])
        load([imageFolder, '/stimuli_sessions.mat'])
        
        % find desired ID of stimulus selection & locations - automatic or
        % user defined
        if sum(strcmp(switch_selection, {'auto', 'automatic'}))
            % find stimulus files with selection numbers
            if strcmp(task_type, 'behavior')
                SelFiles = dir(['./sequences/stim_info_*', 'behavior_*','_sel*','.mat']);
            else
                SelFiles = dir(['./sequences/stim_info_*', '_sel*','_*','.mat']);
            end
            incl = zeros(length(SelFiles),1);
            usedID = [];
            for f = 1:length(SelFiles)
                % is file of the same task type?
                tmp = strfind(SelFiles(f).name,'_sel');
                if ~isempty(tmp)
                    tmp2 = tmp-1+strfind(SelFiles(f).name(tmp:end), '_');
                    usedID = [usedID,str2double(SelFiles(f).name(tmp+4:tmp2(2)-1))];
                end
            end
            freeID = setdiff(1:size(stimuli_sessions,2),usedID);
            selectionID = freeID(1);
        else % user defined
            selectionID = str2double(switch_selection);
        end
        
        % select the desired stimuli and locations
        stimuli_info.session = stimuli_sessions(:,selectionID);
        
        % find positions on the screen and shuffle according to desired ID
        stimuli_info.positions = zeros(size(stimuli_info,1),2);
        for s = 1:nsessions
            [cue_positions,~] = create_circular_array_angleoffset...
                (W/2,H/2,R,angle_offset(s),stimuli_info(stimuli_info.session == s,:));
            locids = mod(stimuli_locations(stimuli_info.session==s,selectionID),nstim);
            locids(locids == 0) = nstim;
            cue_positions_rand = cue_positions(locids,:);
            stimuli_info.positions(stimuli_info.session==s,:) = cue_positions_rand;
        end
        
    end
    
    % save to file
    save(['./sequences/stim_info_',subID,'.mat'], 'stimuli_info')
    save(['./sequences/stim_info_',subID, '_sel',num2str(selectionID),'_' ,now,'.mat'], 'stimuli_info')
else % load existing
    % load stimulus info
    load(['./sequences/stim_info_',subID,'.mat'])
end

% normal session or delayed retest?
if sessionID == 3
    stim_session = ismember(stimuli_info.session, 1);
elseif sessionID == 4
    stim_session = ismember(stimuli_info.session, 2);
else
    stim_session = stimuli_info.session == sessionID;
end

% take only the stimuli for this session
stimuli_info = stimuli_info(stim_session,:);

% take only the locations for this session
cue_positions_rand = stimuli_info.positions;
% order the positions
dum = angle((cue_positions_rand(:,1)-W/2)+1i*(cue_positions_rand(:,2)-H/2));
dum(dum<0) = dum(dum<0)+2*pi;
[~,id] = sort(dum);
cue_positions = cue_positions_rand(id,:);
[~,stim_cue_idx] = ismember(cue_positions_rand,cue_positions, 'rows');

% file names for this session
imgList = stimuli_info.file_name;

%% get sequences of stimuli

if new_participant == 1 || old_session ~= sessionID
    
    if nrep_fam > 0
        % familiarization
        sequence_familiarization = create_sequence_miniblocks(nstim, nrep_fam, ...
            mindiff, nstim*nrep_fam, nrep_fam, 0,0);
        % assign catch question types
        for s = 1:nstim
            sequence_familiarization(sequence_familiarization(:,1)==s,2) = randperm(nrep_fam);
        end
        save(strcat('sequences/m_sequence_familiarization_',subID,'.mat'),'sequence_familiarization');
        save(strcat('sequences/m_sequence_familiarization_',subID,'_',now,'.mat'),'sequence_familiarization');
    end
    
    if ~switch_visual
        
        if nrep_enc > 0
            % encoding
            sequence_encoding = create_sequence_miniblocks(nstim, nrep_enc, ...
                mindiff, ntrials_block_enc, nrep_DDtrials, mindiffDD,1);
            save(strcat('sequences/m_sequence_encoding_',subID,'.mat'),'sequence_encoding');
            save(strcat('sequences/m_sequence_encoding_',subID,'_',now,'.mat'),'sequence_encoding');
        end
        
        if nrep_ret > 0
            % retrieval
            sequence_retrieval = create_sequence_miniblocks(nstim, nrep_ret, ...
                mindiff, ntrials_block_ret, nrep_Ctrials, mindiffC,0);
            % assign catch question types
            for s = 1:nstim
                sequence_retrieval(sequence_retrieval(:,1)==s & sequence_retrieval(:,2) == 1,2) = ...
                    randperm(nrep_Ctrials);
            end
            save(strcat('sequences/m_sequence_retrieval_',subID,'.mat'),'sequence_retrieval');
            save(strcat('sequences/m_sequence_retrieval_',subID,'_',now,'.mat'),'sequence_retrieval');
        end
    end
elseif new_participant == 0
    % load sequence
    if nrep_fam>0
    load(strcat('sequences/m_sequence_familiarization_',subID,'.mat'));
    end
    if ~switch_visual
        if nrep_enc>0
            load(strcat('sequences/m_sequence_encoding_',subID,'.mat'));
        end
        if nrep_ret>0
            load(strcat('sequences/m_sequence_retrieval_',subID,'.mat'));
        end
    end
end

%% prep trigger box

if strcmp(trg_type, 'serial')
    % serial
    trg_handle = IOPort('OpenSerialPort', 'COM3');
elseif strcmp(trg_type, 'ttl')
    % ttl
    trg_handle  = DaqDeviceIndex;
    err = DaqDConfigPort(trg_handle,[],0);
    out_ = DaqDOut(trg_handle,0,0); % reset
elseif sum(strcmp(trg_type, {'lab','labjack'}))
    trg_type = 'labjack';
    lab_init;
    trg_handle = L;
else
    trg_handle = [];
end

%% prep result file

outputfile = fopen([resultsFolder '/resultfile_' num2str(subID) '.txt'],'a');
if outputfile == -1
    error('Cannot open the results file.')
end

% content of log file:
% subID
% sessionID
% block_state: encoding or retrieval
% trial_id: number of trial in the presentation sequence
% trial_type: 0 = normal; 1 = drag-and-drop
% cue_id: location number
% cue_xcoord: (pixels)
% cue_ycoord: (pixels)
% stim_label: e.g. german shepherd
% stim_id: stimulus number
% stim_filename: xxx.jpg
% stim_counter: counter of repetitions for each stimulus
% stim_cat: e.g. dog
% stim_perc1: 1 = drawing; 2 = photograph;
% stim_perc2: 1 = left-facing; 2 = right-facing;
% stim_sem1: 1 = animate; 2 = inanimate
% stim_sem2: 1 = flying; 2 = non-flying;
% RT_encoding (s)
% DaD_resp: 0 = forgotten; 1 = remembered
% RT_DaD_reinst (s)
% DaD_numb_attempts: number of attempt to reach correct location
% RT_DaD_loc (s)
% ret_reinst_resp: 0 = forgotten; 1 = remembered
% RT_ret_reinst
% catch_type: perc1, perc2, sem1, sem2, exemplar
% RT_catch (s)
% catch_resp: 0 = incorrect; 1 = correct; 2 = forgotten; 3 no answer
% onset_session (s)
% onset_familiarization (s)
% onset_encoding (s)
% onset_retrieval (s)
% onset_trial (s)
% onset_cue (s)
% onset_trigger (s)
% onset_stim (s)
% onset_DaD_stim (s)
% onset_DaD_reinst (s)
% onset_DaD_loc (s)
% onset_ret_reinst (s)
% onset_catch (s)
% onset_catch_resp (s)

if new_participant
    fprintf(outputfile, ['subID\tsessionID\tblock_state\ttrial_id\ttrial_type\tcue_id\tcue_xcoord\tcue_ycoord\t',...
        'stim_label\tstim_id\tstim_filename\tstim_counter\tstim_cat\t',...
        'stim_perc1\tstim_perc2\tstim_sem1\tstim_sem2\tRT_encoding\tDaD_resp\tRT_DaD_reinst\tDaD_numb_attempts\tRT_DaD_loc\t',...
        'ret_reinst_resp\tRT_ret_reinst\tcatch_type\tRT_catch\tcatch_resp\t',...
        'onset_session\tonset_familiarization\tonset_encoding\tonset_retrieval\tonset_trial\tonset_cue\tonset_trigger\tonset_stim\t',...
        'onset_DaD_stim\tonset_DaD_reinst\tonset_DaD_loc\tonset_ret_reinst\tonset_catch\tonset_catch_resp\n']);
end


%% INSTRUCTIONS

if sum(strcmp(switch_instructions, {'yes','y','Y','YES'}))
    
    Screen(window1,'FillRect',[1,1,1]);
    Screen('Flip', window1);
    
    if tr_fam==1
        inst_files = dir(['.\instructions\instructions_',task_type,'\instructions_*.jpg']);
        
        ins = 1;
        while ins < length(inst_files)+1
            %         for ins = 1:length(inst_files)
            
            % load instructions
            instr_object = imread(['.\instructions\instructions_',task_type,'\',inst_files(ins).name]);
            instrDisplay = Screen('MakeTexture', window1, instr_object);
            
            instrSize = size(instr_object);
            instrFact = min(H/instrSize(1),W/instrSize(2));
            
            % Calculate image position (center of the screen)
            pos = [(W-instrSize(2)*instrFact)/2 (H-instrSize(1)*instrFact)/2 ...
                (W+instrSize(2)*instrFact)/2 (H+instrSize(1)*instrFact)/2];
            
            Priority(MaxPriority(window1));
            Priority(2);
            
            Screen('DrawTexture', window1, instrDisplay, [], pos);
            Screen('Flip', window1,[],0);
            
            % wait for a key press to continue to the next page
            [~, keyCode, ~] = KbStrokeWait;
            
            if keyCode(KbName(responseKeys{3}))
                ins = ins-1;
                if ins<=0; ins = 1; end
            elseif keyCode(KbName(responseKeys{4}))
                ins = ins+1;
            elseif keyCode(KbName(instrKey{1}))
                ins = ins+1;
            elseif keyCode(KbName(stopKey{1}))
                quitTask; return
            end
        end
    end
    
    Screen(window1,'FillRect',backgroundColor);
    Screen('Flip', window1);
    
end

%% Start Session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SetMouse(0,0,window1)

% time of start of experiment and date
SessionStart = GetSecs - ExpStart;

%% FAMILIARIZATION

if tr_fam<nstim*nrep_fam
    Screen('TextSize', window1, 2*fontsize);
    if switch_visual
        if DPIswitch
            Screen('DrawText', window1, 'VISUAL TASK', W/2-16*fontsize, H/2, textColor);
        else
            DrawFormattedText(window1, 'VISUAL TASK', 'center', H/2, textColor);
        end
    else
        if DPIswitch
            Screen('DrawText', window1, 'FAMILIARIZATION', W/2-16*fontsize, H/2, textColor);
        else
            DrawFormattedText(window1, 'FAMILIARIZATION', 'center', H/2, textColor);
        end
    end
    Screen('Flip', window1);
    WaitSecs(1);
    
    % load instructions
    instr_object = imread(strcat('./instructions/instructions_',task_type,'/reminder_familiarization.jpg'));
    instrDisplay = Screen('MakeTexture', window1, instr_object);
    
    % Calculate image position (center of the screen)
    instrSize = size(instr_object);
    instrFact = min(H/instrSize(1),W/instrSize(2));
    pos = [(W-instrSize(2)*instrFact)/2 (H-instrSize(1)*instrFact)/2 ...
        (W+instrSize(2)*instrFact)/2 (H+instrSize(1)*instrFact)/2];
    Screen('DrawTexture', window1, instrDisplay, [], pos);
    
    %     Screen('TextSize', window1, 20);
    %     Screen('DrawText',window1, 'You will see an object on the screen, followed by two possible descriptions.',instroffset,H/5, black);
    %     Screen('DrawText',window1, 'Please choose the correct description by pressing LEFT or RIGHT.',instroffset,H/5+1*H/5, black);
    %     Screen('DrawText',window1, 'You will then see whether your answer was correct or incorrect.',instroffset,H/5+2*H/5, black);
    %     Screen('DrawText',window1, 'Press any key to continue.',instroffset,H/5+3*H/5, black);
    
    Screen('Flip',window1);
    KbStrokeWait;
end

% send start trigger
for i = 1:3
    sendTrigger(trg_type, trg_handle);
end


block_state = 'familiarization';
FamiliarizationStart = GetSecs - ExpStart;

% loop over trials
for t = tr_fam:nrep_fam*nstim
    
    % reset time and RT variables
    [RT_familiarization,catch_type, catch_resp, ...
        TrialStart, TriggerStart, CatchStart, StimStart, FamRespStart ] = deal([]);
    
    KbReleaseWait();
    
    % PREP stimuli
    file_object_1 = imgList{(sequence_familiarization(t,1))};
    img_object_1 = imread(fullfile(imageFolder ,file_object_1));
    objectDisplay_1 = Screen('MakeTexture', window1, img_object_1);
    
    % Calculate image position (center of the screen) OBJECTS
    pos_o = [(W-imageSize_o(2))/2 (H-imageSize_o(1))/2 (W+imageSize_o(2))/2 (H+imageSize_o(1))/2];
    
    % PREP catch quesion
    % get catch_type;
    dum = mod(sequence_familiarization(t,2),length(catch_types));
    if dum == 0; dum = length(catch_types); end
    catch_type = catch_types{dum};
    
    % find labels
    correct_label = eval(['stimuli_info.label_',catch_type,'(sequence_familiarization(t,1))']);
    if strcmp(catch_type,'exemplar') && length(unique(stimuli_info.cat_exemplar))<nstim
        catids = find(stimuli_info.cat_exemplar == stimuli_info.cat_exemplar(sequence_familiarization(t,1)));
        other_labels = stimuli_info.label_exemplar(setdiff(catids,sequence_familiarization(t,1)));
    else
        other_labels = setdiff(...
            unique(eval(['stimuli_info.label_',catch_type,'(setdiff(1:nstim,sequence_familiarization(t,1)),:)'])),...
            correct_label);
    end
    labels(1) = correct_label;
    labels(2) = other_labels(randperm(length(other_labels),1)); % in case there is more than one lures, pick one at random
    
    % randomize location of labels
    xoffset = [length(labels{1}),length(labels{2})]*fontsize/2;
    rand_positions = catch_positions(randperm(size(catch_positions,1)),:);
    
    % --------- Start trial
    TrialStart = GetSecs - ExpStart;
    
    % Screen priority
    Priority(MaxPriority(window1));
    Priority(2);
    
    
    % ---------- Show fixation cross + arena
    Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
    drawCross(window1,W,H,barColor);
    tFixation = Screen('Flip', window1,[],1); % use 1 here to allow addition of the location!
    Screen('Flip', window1, tFixation + (fixationDuration + rand(1))- slack,1); %Adding between 0 and 1 second randomly to the minimum for the fixation duration
    
    % if behavioral task, show question first
    if switch_behavior || switch_visual
        Screen('TextSize', window1, fontsize);
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        drawCross(window1,W,H,barColor);
        if DPIswitch
            Screen('DrawText', window1, labels{1},(rand_positions(1,1)) - xoffset(1) +labeloffset, catch_y, textColor);%CORRECT
            Screen('DrawText', window1, labels{2},(rand_positions(2,1)) - xoffset(2) +labeloffset, catch_y, textColor);%LURE
        else
            DrawFormattedText_mod(window1, labels{1},'center', catch_y,textColor,[],[],[],[],[],[],(rand_positions(1,1)));%CORRECT
            DrawFormattedText_mod(window1, labels{2},'center', catch_y, textColor,[],[],[],[],[],[],(rand_positions(2,1)));%LURE
        end
        
        Screen('Flip', window1,[],1);
        WaitSecs(catchTimeout);
    end
    
    % -------- Send Trigger
    TriggerStart = GetSecs - ExpStart;
    sendTrigger(trg_type, trg_handle);
    
    % ------ show stimulus
    %     Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
    %     Screen('Flip', window1,0.5,1);
    Screen('DrawTexture', window1, objectDisplay_1, [], pos_o);
    tdum = Screen('Flip', window1,[],1);
    FamStimStart = tdum - ExpStart;
    if switch_behavior || switch_visual
        CatchStart = FamStimStart;
    end
    
    % if standard task, show question now
    if switch_standard
        WaitSecs(famTimeout);
        
        % ------ show question after object
        Screen('Flip',window1,slack,0);
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        Screen('TextSize', window1, fontsize);
        if DPIswitch
            Screen('DrawText', window1, labels{1},(rand_positions(1,1)) - xoffset(1) +labeloffset, catch_y, textColor);%CORRECT
            Screen('DrawText', window1, labels{2},(rand_positions(2,1)) - xoffset(2) +labeloffset, catch_y, textColor);%LURE
        else
            DrawFormattedText_mod(window1, labels{1},'center', catch_y,textColor,[],[],[],[],[],[],(rand_positions(1,1)));%CORRECT
            DrawFormattedText_mod(window1, labels{2},'center', catch_y, textColor,[],[],[],[],[],[],(rand_positions(2,1)));%LURE
        end
        tdum = Screen('Flip', window1,[],1);
        CatchStart = tdum - ExpStart;
    end
    
    % ------ wait for response
    % wait for response
    waitForResponse = true;
    while waitForResponse
        [keyIsDown,respTime,keyCode] = KbCheck;
        pressedKeys = find(keyCode);
        
        if ~isempty(pressedKeys) && sum(ismember(KbName(responseKeys(3:end)),pressedKeys(1)))
            catch_resp_key = responseKeys{find(ismember(KbName(responseKeys(3:end)),pressedKeys(1)))+2};
            CatchRespStart = respTime - ExpStart;
            RT_familiarization = CatchRespStart - CatchStart;
            waitForResponse = false;
        elseif keyCode(KbName(stopKey{1}))
            quitTask; return
        end
    end
    
    % check answer
    catch_resp = 2; % no answer
    if switch_visual
        catch_resp = 2; % no answer
        Screen('TextSize', window1, fontsize);
        if (rand_positions(1,1)==catch_positions(1,1) && catch_resp_key(1)==responseKeys{3}(1)) || (rand_positions(1,1)==catch_positions(2,1) && catch_resp_key(1)==responseKeys{4}(1))
            catch_resp = 1; %correct
        elseif (rand_positions(1,1)==catch_positions(1,1) && catch_resp_key(1)==responseKeys{4}(1)) || (rand_positions(1,1)==catch_positions(2,1) && catch_resp_key(1)==responseKeys{3}(1))%'R'
            catch_resp = 0; %incorrect
        end
        Screen('Flip',window1,0,0);
    else
        Screen('TextSize', window1, fontsize);
        if rand_positions(1,1)==catch_positions(1,1) && catch_resp_key(1)==responseKeys{3}(1) %'L'
            %             correct_responses_catch = correct_responses_catch + 1;
            catch_resp = 1; %correct
            if DPIswitch
                Screen('DrawText', window1, labels{1},(rand_positions(1,1)) - xoffset(1)+labeloffset, catch_y, barColor_Correct);%CORRECT
            else
                DrawFormattedText_mod(window1, labels{1},'center', catch_y, barColor_Correct,[],[],[],[],[],[],(rand_positions(1,1)));%CORRECT
            end
        elseif rand_positions(1,1)==catch_positions(1,1) && catch_resp_key(1)==responseKeys{4}(1) %'R'
            catch_resp = 0; %incorrect
            if DPIswitch
                Screen('DrawText', window1, labels{2},(rand_positions(2,1)) - xoffset(2)+labeloffset, catch_y, barColor_Incorrect);
            else
                DrawFormattedText_mod(window1, labels{2},'center', catch_y,barColor_Incorrect,[],[],[],[],[],[],(rand_positions(2,1)));
            end
        elseif rand_positions(1,1)==catch_positions(2,1) && catch_resp_key(1)==responseKeys{4}(1) %'R'
            %             correct_responses_catch = correct_responses_catch + 1;
            catch_resp = 1;
            if DPIswitch
                Screen('DrawText', window1, labels{1},(rand_positions(1,1)) - xoffset(1)+labeloffset, catch_y, barColor_Correct);%CORRECT
            else
                DrawFormattedText_mod(window1, labels{1},'center', catch_y,barColor_Correct,[],[],[],[],[],[],(rand_positions(1,1)));
            end
        elseif rand_positions(1,1)==catch_positions(2,1) && catch_resp_key(1)==responseKeys{3}(1) %'L'
            catch_resp = 0;
            if DPIswitch
                Screen('DrawText', window1, labels{2},(rand_positions(2,1)) - xoffset(2)+labeloffset, catch_y, barColor_Incorrect);
            else
                DrawFormattedText_mod(window1, labels{2},'center', catch_y,barColor_Incorrect,[],[],[],[],[],[],(rand_positions(2,1)));
            end
        end
        Screen('Flip',window1,0,0);
        WaitSecs(feedbacktimeout);
    end
    
    % ---------- print everything to log file
    stimID = sequence_familiarization(t,1);
    
    fprintf(outputfile, ['%s\t %d\t %s\t %d\t %d\t %d\t %f\t %f\t',...
        '%s\t %d\t %s\t %d\t %d\t ',...
        '%d\t %d\t %d\t %d\t %f\t %d\t %f\t %d\t %f\t',...
        '%d\t %f\t %s\t %f\t %d\t',...
        '%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t',...
        '%f\t %f\t %f\t %f\t %f\t %f\n'] ,...
        subID, sessionID, block_state, t, sequence_familiarization(t,2), stim_cue_idx(stimID), stimuli_info.positions(stimID,1), stimuli_info.positions(stimID,2),...
        stimuli_info.item_name{stimID}, stimuli_info.item_ID(stimID,:), imgList{stimID}, sum(sequence_familiarization(1:t,1)==stimID), stimuli_info.cat_exemplar(stimID),...
        stimuli_info.cat_perc_1(stimID), stimuli_info.cat_perc_2(stimID), stimuli_info.cat_sem_1(stimID),stimuli_info.cat_sem_2(stimID),...
        [], [], [], [], [],...
        [],[],catch_type,RT_familiarization,catch_resp,...
        SessionStart,FamiliarizationStart,[],[],TrialStart,[],TriggerStart,FamStimStart,...
        [],[],[], [],CatchStart,CatchRespStart);
end

if switch_visual
    if DPIswitch
        Screen('DrawText',window1, 'Well done!',W/2-2*50,H/4, textColor);
    else
        DrawFormattedText(window1, 'Well done!','center',H/4, textColor);
    end
    Screen('Flip',window1,1,[]);
    
    KbStrokeWait;
    KbReleaseWait();
    RestrictKeysForKbCheck([]);
    Screen(window1,'Close');
    fclose('all');
    close all
    sca;
    return
end

%% ENCODING

if tr_enc<nstim*nrep_enc
    Screen('TextSize', window1, 2*fontsize);
    if DPIswitch
        Screen('DrawText', window1, 'ENCODING', W/2-9*fontsize, H/2, textColor);
    else
        DrawFormattedText(window1, 'ENCODING', 'center', H/2, textColor);
    end
    % DrawFormattedText(window1, 'ENCODING','center',H/2, black);
    Screen('Flip',window1);
    WaitSecs(1);
    
    % load instructions
    instr_object = imread(strcat('./instructions/instructions_',task_type,'/reminder_encoding.jpg'));
    instrDisplay = Screen('MakeTexture', window1, instr_object);
    
    % Calculate image position (center of the screen)
    instrSize = size(instr_object);
    instrFact = min(H/instrSize(1),W/instrSize(2));
    pos = [(W-instrSize(2)*instrFact)/2 (H-instrSize(1)*instrFact)/2 ...
        (W+instrSize(2)*instrFact)/2 (H+instrSize(1)*instrFact)/2];
    Screen('DrawTexture', window1, instrDisplay, [], pos);
    
    %     Screen('TextSize', window1, 20);
    %     Screen('DrawText',window1, 'One location around the circle will be highlighted, then an object is shown.',instroffset,H/6, black);
    %     Screen('DrawText',window1, 'When you remember the location of the object, press UP.',instroffset,H/6+1*H/8, black);
    %     Screen('DrawText',window1, 'Every now and then, you are shown an object and asked to indicate their location.',instroffset,H/6+2*H/8, black);
    %     Screen('DrawText',window1, 'Use LEFT and RIGHT to move around the circle and use UP to confirm a location.',instroffset,H/6+3*H/8, black);
    %     Screen('DrawText',window1, 'You are shown whether your answer was correct or incorrect.',instroffset,H/6+4*H/8, black);
    %     Screen('DrawText',window1, 'Press any key to continue.',instroffset,H/6+5*H/8, black);
    
    Screen('Flip',window1);
    KbStrokeWait;
end

% send start trigger
for i = 1:4
    sendTrigger(trg_type, trg_handle);
end

block_state = 'encoding';
EncodingStart = GetSecs - ExpStart;

% loop over trials
for t = tr_enc:nrep_enc*nstim
    
    % reset time and RT variables
    [RT_encoding, DaD_resp, RT_DaD_reinst, DaD_numb_attempts, RT_DaD_loc,...
        TrialStart,CueStart,TriggerStart,StimStart,...
        DaDStimStart,DaDReinstStart,DaDLocRespStart] = deal([]);
    
    KbReleaseWait();
    
    % check for mini-break
    if t>tr_enc && mod(t,ntrials_block_enc) == 1
        
        % Fill up a bar during the break
        Screen('TextSize', window1, fontsize);
        centeredRect = CenterRectOnPointd([0 0 barWidth+20 barHeight+20], xCenter,(H/1.5));
        growingsize = barWidth/(hz*minBreakTime);
        
        for fill = 1:hz*minBreakTime % fill up the bar
            if DPIswitch
                Screen('DrawText', window1, 'Time for a mini break!', W/2-22*fontsize/2, H/3, textColor);
            else
                DrawFormattedText(window1, 'Time for a mini break!', 'center', H/3, textColor);
            end
            Screen('FillRect', window1, [255 255 255], centeredRect);
            Screen('FillRect', window1, [200 200 200], [xCenter-barWidth/2 (H/1.5)-barHeight/2 xCenter-barWidth/2+growingsize*fill (H/1.5)+barHeight/2]);%centeredRectfill);
            Screen('Flip', window1,1);
        end
        if DPIswitch
            Screen('DrawText', window1, 'When you''re ready, press an arrow key to continue.', W/2-53*fontsize/2, H/3, textColor);
        else
            DrawFormattedText(window1, 'When you''re ready, press an arrow key to continue.','center',H/3, textColor);
        end
        Screen('Flip',window1);
        KbStrokeWait;
    end
    
    KbReleaseWait();
    
    % prep stimuli
    file_object_1 = imgList{(sequence_encoding(t,1))};
    img_object_1 = imread(fullfile(imageFolder ,file_object_1));
    objectDisplay_1 = Screen('MakeTexture', window1, img_object_1);
    
    % Calculate image position (center of the screen) OBJECTS
    pos_o = [(W-imageSize_o(2))/2 (H-imageSize_o(1))/2 (W+imageSize_o(2))/2 (H+imageSize_o(1))/2];
    
    TrialStart = GetSecs - ExpStart;
    
    % test for normal presentation or position test trial
    if sequence_encoding(t,2) == 1
        % do a drag-and-drop trial!!
        
        % Screen priority
        Priority(MaxPriority(window1));
        Priority(2);
        
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        Screen('Flip', window1,[],1);
        
        % instruction text (self paced)
        Screen('TextSize', window1, fontsize);
        if DPIswitch
            Screen('DrawText', window1, 'Select the location of the next object. Press a key to continue...', W/2-67*fontsize/2, (H-R)/4, textColor);
        else
            DrawFormattedText(window1, 'Select the location of the next object. Press a key to continue...','center',(H-R)/4, textColor);
        end
        Screen('Flip',window1)
        WaitSecs(0.2);
        [~, keyCode, ~] = KbStrokeWait;
        
        % -------- Send Trigger
        TriggerStart = GetSecs - ExpStart;
        sendTrigger(trg_type, trg_handle)
        
        % show object
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        %         Screen('Flip', window1,0.5,1);
        Screen('DrawTexture', window1, objectDisplay_1, [], pos_o);
        tdum = Screen('Flip', window1,[],0); % Start of drag & drop trial
        DaDStimStart = tdum - ExpStart;
        
        % wait for reinstatement button press
        waitForResponse = true;
        while waitForResponse
            [~, DaDReinstStart, keyCode] = KbCheck;
            % Check for response keys
            if keyCode(KbName(responseKeys{1})) %keyCode(KbName('UpArrow')) % remembered
                DaD_resp = 1;
                RT_DaD_reinst = DaDReinstStart -ExpStart - DaDStimStart;
                waitForResponse = false;
            elseif keyCode(KbName(responseKeys{2}))% keyCode(KbName('DownArrow')) % forgotten
                DaD_resp = 0;
                RT_DaD_reinst = DaDReinstStart -ExpStart - DaDStimStart;
                waitForResponse = false;
            elseif keyCode(KbName(stopKey{1}))
                quitTask; return
            end
        end
        
        KbReleaseWait();
        clear keyCode
        
        % show locations
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        Screen('DrawDots', window1, stimuli_info.positions', dotsize, black, [], []);
        % Screen('Flip',window1,[],1);
        
        % start at random location
        correct_loc = find(sum([(cue_positions(:,1) - stimuli_info.positions(sequence_encoding(t,1),1))==0,...
            (cue_positions(:,2) - stimuli_info.positions(sequence_encoding(t,1),2))==0],2) == 2);
        start_locs = cue_positions(setdiff(1:nstim,correct_loc),:);
        rand_loc = randi(nstim-1,1);
        Screen('FillOval', window1, barColorDD,[(start_locs(rand_loc,1)-2*dotsize) (start_locs(rand_loc,2)-2*dotsize) (start_locs(rand_loc,1)+2*dotsize) (start_locs(rand_loc,2)+2*dotsize)]);
        tdum = Screen('Flip',window1,[],0);
        DaDLocStart = tdum - ExpStart;
        
        % move cursor
        contDD = true;
        loc_ids = setdiff(1:nstim,correct_loc);
        cur_loc = loc_ids(rand_loc); % keep track of cursor location using cur_loc
        DaD_numb_attempts = 0;
        while contDD
            [keyIsDown,respTime,keyCode] = KbCheck;
            
            if keyCode(KbName(responseKeys{3})) % move counter-clockwise
                if cur_loc == 1
                    new_loc = nstim;
                else
                    new_loc = cur_loc-1;
                end
                Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
                Screen('DrawDots', window1, cue_positions', dotsize, black, [], []);
                Screen('FillOval', window1, barColorDD,[(cue_positions(new_loc,1)-2*dotsize) (cue_positions(new_loc,2)-2*dotsize) (cue_positions(new_loc,1)+2*dotsize) (cue_positions(new_loc,2)+2*dotsize)]);
                Screen('Flip',window1,[],0);
                cur_loc = new_loc;
                WaitSecs(0.2)
            elseif keyCode(KbName(responseKeys{4})) % move clockwise
                if cur_loc == nstim
                    new_loc = 1;
                else
                    new_loc = cur_loc+1;
                end
                Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
                Screen('DrawDots', window1, stimuli_info.positions', dotsize, black, [], []);
                Screen('FillOval', window1, barColorDD,[(cue_positions(new_loc,1)-2*dotsize) (cue_positions(new_loc,2)-2*dotsize) (cue_positions(new_loc,1)+2*dotsize) (cue_positions(new_loc,2)+2*dotsize)]);
                Screen('Flip',window1,[],0);
                cur_loc = new_loc;
                WaitSecs(0.2)
            elseif keyCode(KbName(responseKeys{1}))% keyCode(KbName('UpArrow'))
                DaD_numb_attempts = DaD_numb_attempts + 1;
                % check answer and give feedback
                if sum((cue_positions(cur_loc,:) - stimuli_info.positions(sequence_encoding(t,1),:))==0) == 2
                    % correct
                    Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
                    Screen('DrawDots', window1, stimuli_info.positions', dotsize, black, [], []);
                    Screen('FillOval', window1, barColor_Correct,[(cue_positions(cur_loc,1)-2*dotsize) (cue_positions(cur_loc,2)-2*dotsize) (cue_positions(cur_loc,1)+2*dotsize) (cue_positions(cur_loc,2)+2*dotsize)]);
                    feedbackStart = Screen('Flip',window1,[],1);
                    DaDLocRespStart = respTime - ExpStart;
                    RT_DaD_loc = DaDLocRespStart - DaDLocStart;
                    Screen('Flip',window1,feedbackStart+feedbacktimeout,0);
                    contDD = false;
                else
                    % incorrect
                    Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
                    Screen('DrawDots', window1, stimuli_info.positions', dotsize, black, [], []);
                    Screen('FillOval', window1, barColor_Incorrect,[(cue_positions(cur_loc,1)-2*dotsize) (cue_positions(cur_loc,2)-2*dotsize) (cue_positions(cur_loc,1)+2*dotsize) (cue_positions(cur_loc,2)+2*dotsize)]);
                    Screen('Flip',window1,[],0);
                    WaitSecs(0.2)
                end
            elseif keyCode(KbName(stopKey{1}))
                quitTask; return
            end
        end
        
        % instructions
        %         Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        %         Screen('Flip', window1,[],1);
        %         Screen('TextSize', window1, 20);
        %         DrawFormattedText(window1, 'Well done!','center',H/5, black);
        %         DrawFormattedText(window1, 'We will continue with a normal trial.','center',4*H/5, black);
        %         Screen('Flip', window1,[],1);
        %         KbStrokeWait;
        
    else
        % normal presentation
        
        % Screen priority
        Priority(MaxPriority(window1));
        Priority(2);
        
        % ---------- Show fixation cross + arena
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        Screen('DrawDots', window1, stimuli_info.positions', dotsize, black, [], []);
        drawCross(window1,W,H,barColor);
        tFixation = Screen('Flip', window1,[],1); % use 1 here to allow addition of the location!
        Screen('Flip', window1, tFixation + (fixationDuration + rand(1))- slack,1); %Adding between 0 and 1 second randomly to the minimum for the fixation duration
        
        % ----------- Show CUE
        drawCross(window1,W,H,barColor);
        Screen('FillOval', window1, barColorCue,[ (stimuli_info.positions(sequence_encoding(t,1),1)-2*dotsize) (stimuli_info.positions(sequence_encoding(t,1),2)-2*dotsize) (stimuli_info.positions(sequence_encoding(t,1),1)+2*dotsize) (stimuli_info.positions(sequence_encoding(t,1),2)+2*dotsize)]);
        
        % start of cue period
        tdum = Screen('Flip', window1); % Start of trial
        CueStart = tdum - ExpStart;
        
%         while GetSecs - (CueStart+ExpStart) < cueTimeout+rand(1)
%             [keyIsDown,secs,keyCode] = KbCheck;
%             
%             % check for quit key press
%             if keyCode(KbName(stopKey{1}))
%                 quitTask; return
%             end
%         end
        WaitSecs(cueTimeout+rand(1))
        
        % ------- Show OBJECTs
        % Screen priority
        Priority(MaxPriority(window1));
        Priority(2);
        
        % -------- Show fixation cross
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        drawCross(window1,W,H,barColor);
        tFixation = Screen('Flip', window1,[],1);
        Screen('Flip', window1, tFixation + (fixationDuration)- slack,0);
        
        % -------- Send Trigger
        TriggerStart = GetSecs - ExpStart;
        sendTrigger(trg_type, trg_handle)
        
        % -------- Show Object
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        Screen('DrawTexture', window1, objectDisplay_1, [], pos_o);
        tdum = Screen('Flip', window1); % Start of trial
        StimStart = tdum - ExpStart;
        
        waitForResp = true;
        while waitForResp || GetSecs - (StimStart+ExpStart) < objectTimeout
            [keyIsDown,respTime,keyCode] = KbCheck;
            
            % ESC key quits the experiment
            if keyCode(KbName(stopKey{1}))
                quitTask; return
            end
            
            % Check for response keys
            if sum(keyCode)
                RT_encoding = respTime - ExpStart - StimStart;
                waitForResp = false;
            end
        end
        
    end
    
    % --------- inter trial interval
    Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
    ITIStart = Screen('Flip', window1);
    
    while GetSecs - ITIStart < postobjectTimeout
        [keyIsDown,secs,keyCode] = KbCheck;
        
        % ESC key quits the experiment
        if keyCode(KbName(stopKey{1}))
            quitTask; return
        end
    end
    
    % ---------- print everything to log file
    stimID = sequence_encoding(t,1);
    
    fprintf(outputfile, ['%s\t %d\t %s\t %d\t %d\t %d\t %f\t %f\t',...
        '%s\t %d\t %s\t %d\t %d\t ',...
        '%d\t %d\t %d\t %d\t %f\t %d\t %f\t %d\t %f\t',...
        '%d\t %f\t %s\t %f\t %d\t',...
        '%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t',...
        '%f\t %f\t %f\t %f\t %f\t %f\n'] ,...
        subID, sessionID, block_state, t, sequence_encoding(t,2), stim_cue_idx(stimID), stimuli_info.positions(stimID,1), stimuli_info.positions(stimID,2),...
        stimuli_info.item_name{stimID}, stimuli_info.item_ID(stimID,:), file_object_1, sum(sequence_encoding(1:t,1)==stimID), stimuli_info.cat_exemplar(stimID),...
        stimuli_info.cat_perc_1(stimID), stimuli_info.cat_perc_2(stimID), stimuli_info.cat_sem_1(stimID),stimuli_info.cat_sem_2(stimID),...
        RT_encoding, DaD_resp, RT_DaD_reinst, DaD_numb_attempts, RT_DaD_loc,...
        [],[],[],[],[],...
        SessionStart,[], EncodingStart,[],TrialStart,CueStart,TriggerStart,StimStart,...
        DaDStimStart,DaDReinstStart,DaDLocRespStart, [],[],[]);
end

%% DISTRACTOR

if tr_ret == 1 && distr_duration > 0
    Screen('TextSize', window1, 2*fontsize);
    if DPIswitch
        Screen('DrawText', window1, 'NUMBER TASK', W/2-18*fontsize, H/2, textColor);
    else
        DrawFormattedText(window1, 'NUMBER TASK','center',H/2, textColor);
    end
    Screen('Flip',window1);
    WaitSecs(1);
    
    % INSTRUCTION DISTRACTOR
    % load instructions
    instr_object = imread(strcat('./instructions/instructions_',task_type,'/reminder_distractor.jpg'));
    instrDisplay = Screen('MakeTexture', window1, instr_object);
    
    % Calculate image position (center of the screen)
    instrSize = size(instr_object);
    instrFact = min(H/instrSize(1),W/instrSize(2));
    pos = [(W-instrSize(2)*instrFact)/2 (H-instrSize(1)*instrFact)/2 ...
        (W+instrSize(2)*instrFact)/2 (H+instrSize(1)*instrFact)/2];
    Screen('DrawTexture', window1, instrDisplay, [], pos);
    
    %     Screen('TextSize',window1,18);
    %     Screen('DrawText',window1, 'Is the number displayed ODD or EVEN?', instroffset, (H/4), [0 0 0]);
    %     distractor_message = strcat('Try to complete as many trials as you can in (',num2str(distr_duration),') seconds');
    %     Screen('DrawText',window1, distractor_message, instroffset,H/4+H/4, [0 0 0]);
    %     Screen('DrawText',window1, 'Use <- for ODD and -> for EVEN',instroffset, H/4+2*H/4, [0 0 0]);
    %     Screen('DrawText',window1, 'Press any key to begin',instroffset, H/4+3*H/4, [0 0 0]);
    
    Screen('Flip',window1);
    Screen('TextSize',window1,2*fontsize);
    
    % Wait for subject to press a key
    WaitSecs(0.2);
    [secs,keyCode,~] = KbWait;
    if keyCode(KbName(stopKey{1}))
        quitTask; return
    end
    clear keyCode
    
    correct_responses_distractor = 0;
    number_trials_distractor = 0;
    
    KbReleaseWait();
    
    startTime_distractor_task = Screen('Flip', window1);
    while GetSecs - startTime_distractor_task < distr_duration
        %Load numbers for the distractor task
        number = randi([10,99]);
        correct_response_distractor = rem((number),2); %0 = EVEN, 1=ODD
        
        %Convert to string the number for the distractor trial
        mes_distractor = strcat(num2str(number));
        
        % Screen priority
        Priority(MaxPriority(window1));
        Priority(2);
        
        % Show fixation cross
        drawCross(window1,W,H);
        tFixation = Screen('Flip', window1);
        
        % Show fixation cross
        Screen(window1, 'FillRect', backgroundColor);
        Screen('Flip', window1, tFixation + fixationDuration - slack,0);
        
        %Showing the numbers
        Screen('TextSize',window1,2*fontsize);
        if DPIswitch
            Screen('DrawText',window1, mes_distractor, W/2-60, H/2-60, [0 0 0]);
        else
            DrawFormattedText(window1, mes_distractor, 'center', 'center', [0 0 0]);
        end
        startTime = Screen('Flip', window1);
        
        % Get keypress response for distractor
        rt_text = 0;
        resp_distractor = 0;
        waitForResponse = true;
        while waitForResponse
            [keyIsDown,secs,keyCode] = KbCheck;
            respTime = GetSecs;
            
            if keyCode(KbName(responseKeys{4}))
                resp_distractor = 'O'; % odd
                rt_text = respTime - startTime;
                waitForResponse = false;
            elseif keyCode(KbName(responseKeys{3}))
                resp_distractor = 'E'; % even
                rt_text = respTime - startTime;
                waitForResponse = false;
            elseif keyCode(KbName(stopKey{1}))
                quitTask; return
            end
        end
        clear keyCode
        
        %Counting the number of responses
        number_trials_distractor = number_trials_distractor + 1;
        
        %Obtaining the accuracy during distractor
        if resp_distractor(1)=='O' && correct_response_distractor ~= 0
            correct_responses_distractor = correct_responses_distractor + 1;
        elseif resp_distractor(1)=='E' && correct_response_distractor == 0
            correct_responses_distractor = correct_responses_distractor + 1;
        end
    end
    
    %%% GIVING CHEER UP MESSAGE OR FEEDBACK FOR DISTRACTOR
    cheer_up_message = {'WELL DONE!','GREAT!','GOOD JOB!'};
    
    Screen('TextSize',window1,fontsize);
    if DPIswitch
        Screen('DrawText',window1, cheer_up_message{round((2)*rand+1)}, W/2-12*fontsize, H/3, [0 0 0]);
        Screen('DrawText',window1, 'Press any key to continue', W/2-6*50, (H/1.5), [0 0 0]);
    else
        DrawFormattedText(window1, cheer_up_message{round((2)*rand+1)}, 'center', H/3, [0 0 0]);
        DrawFormattedText(window1, 'Press any key to continue', 'center', (H/1.5), [0 0 0]);
    end
    Screen('Flip',window1);
    
    % Wait for subject to press a key
    KbStrokeWait;
    
end

%% RETRIEVAL

Screen('TextSize', window1, 2*fontsize);
if DPIswitch
    Screen('DrawText',window1, 'RETRIEVAL',W/2-10*fontsize,H/2, textColor);
else
    DrawFormattedText(window1, 'RETRIEVAL','center',H/2, textColor);
end
Screen('Flip',window1);
WaitSecs(1);

% load instructions
if sessionID <3
    instr_object = imread(strcat('./instructions/instructions_',task_type,'/reminder_retrieval.jpg'));
else
    instr_object = imread(strcat('./instructions/instructions_',task_type,'/reminder_retest.jpg'));
end
instrDisplay = Screen('MakeTexture', window1, instr_object);

% Calculate image position (center of the screen)
instrSize = size(instr_object);
instrFact = min(H/instrSize(1),W/instrSize(2));
pos = [(W-instrSize(2)*instrFact)/2 (H-instrSize(1)*instrFact)/2 ...
    (W+instrSize(2)*instrFact)/2 (H+instrSize(1)*instrFact)/2];
Screen('DrawTexture', window1, instrDisplay, [], pos);

% Screen('TextSize', window1, 20);
% Screen('DrawText',window1, 'One location around the circle will be highlighted.',instroffset,H/6, black);
% Screen('DrawText',window1, 'When you remember the corresponding object, press the UP arrow.',instroffset,H/6+1*H/6, black);
% Screen('DrawText',window1, 'Then think about the object for 3 seconds. Try to remember the details.',instroffset,H/6+2*H/6, black);
% Screen('DrawText',window1, 'Every now and then, you will be shown two possible descriptions of the object.',instroffset,H/6+3*H/6, black);
% Screen('DrawText',window1, 'Choose the correct description by pressing LEFT or RIGHT.',instroffset,H/6+4*H/6, black);
% Screen('DrawText',window1, 'Press any key to continue.',instroffset,H/6+5*H/6, black);

Screen('Flip',window1);
KbStrokeWait;

% send start trigger
for i = 1:5
    sendTrigger(trg_type, trg_handle);
end

block_state = 'retrieval';
RetrievalStart = GetSecs - ExpStart;

for t = tr_ret:nrep_ret*nstim
    
    % reset time and RT variables
    [ret_reinst_resp,RT_ret_reinst,catch_type,RT_catch,catch_resp,...
        TrialStart,CueStart,TriggerStart,RetReinstStart,CatchStart,...
        CatchRespStart] = deal([]);
    
    KbReleaseWait();
    
    % check for mini-break
    if t>tr_ret && mod(t,ntrials_block_ret) == 1
        % Fill up a bar during the break
        Screen('TextSize', window1, fontsize);
        centeredRect = CenterRectOnPointd([0 0 barWidth+20 barHeight+20], xCenter,(H/1.5));
        growingsize = barWidth/(hz*minBreakTime);
        for fill = 1:hz*minBreakTime % fill up the bar
            if DPIswitch
                Screen('DrawText',window1, 'Time for a mini break!',W/2-22*fontsize/2,H/3, textColor);
            else
                DrawFormattedText(window1, 'Time for a mini break!','center',H/3, textColor);
            end
            Screen('FillRect', window1, [255 255 255], centeredRect);
            Screen('FillRect', window1, [200 200 200], [xCenter-barWidth/2 (H/1.5)-barHeight/2 xCenter-barWidth/2+growingsize*fill (H/1.5)+barHeight/2]);%centeredRectfill);
            Screen('Flip', window1,1);
        end
        if DPIswitch
            Screen('DrawText',window1, 'When you''re ready, press an arrow key to continue.',W/2-53*fontsize/2,H/3, textColor);
        else
            DrawFormattedText(window1, 'When you''re ready, press an arrow key to continue.','center',H/3, textColor);
        end
        Screen('Flip',window1);
        KbStrokeWait;
        KbReleaseWait();
    end
    
    % PREP catch question
    if sequence_retrieval(t,2) > 0
        % get catch_type;
        dum = mod(sequence_retrieval(t,2),length(catch_types));
        if dum == 0; dum = length(catch_types); end
        catch_type = catch_types{dum};
        
        % find labels
        correct_label = eval(['stimuli_info.label_',catch_type,'(sequence_retrieval(t,1))']);
        if strcmp(catch_type,'exemplar') && length(unique(stimuli_info.cat_exemplar))<nstim
            catids = find(stimuli_info.cat_exemplar == stimuli_info.cat_exemplar(sequence_retrieval(t,1)));
            other_labels = stimuli_info.label_exemplar(setdiff(catids,sequence_retrieval(t,1)));
        else
            other_labels = setdiff(...
                unique(eval(['stimuli_info.label_',catch_type,'(setdiff(1:nstim,sequence_retrieval(t,1)),:)'])),...
                correct_label);
        end
        labels(1) = correct_label;
        labels(2) = other_labels(randperm(length(other_labels),1)); % in case there is more than one lures, pick one at random
        
        % randomize location of labels
        xoffset = [length(labels{1}),length(labels{2})]*fontsize/2;
        rand_positions = catch_positions(randperm(size(catch_positions,1)),:);
    end
    
    % ------ Start Trial
    TrialStart = GetSecs - ExpStart;
    
    % Screen priority
    Priority(MaxPriority(window1));
    Priority(2);
    
    if sequence_retrieval(t,2) > 0 && switch_behavior
        % ask catch question before reinstatement
        
        % show labels
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        drawCross(window1,W,H,barColor);
        Screen('TextSize', window1, fontsize);
        if DPIswitch
            Screen('DrawText', window1, labels{1},(rand_positions(1,1)) - xoffset(1)+labeloffset, catch_y, textColor);%CORRECT
            Screen('DrawText',window1, '(forgotten)', W/2-100, catch_y+fontsize, [0 0 0]);
            Screen('DrawText', window1, labels{2},(rand_positions(2,1)) - xoffset(2)+labeloffset, catch_y, textColor);%LURE
        else
            DrawFormattedText_mod(window1, labels{1},'center', catch_y,textColor,[],[],[],[],[],[],(rand_positions(1,1)));
            DrawFormattedText_mod(window1, labels{2},'center', catch_y,textColor,[],[],[],[],[],[],(rand_positions(2,1)));
            DrawFormattedText(window1, '(forgotten)','center', catch_y+2*fontsize,textColor,[],[],[],[],[],[]);
        end
        Screen('Flip', window1,0,1);
        WaitSecs(catchTimeout);
    end
    
    % ---------- Show fixation cross + arena
    Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
    Screen('DrawDots', window1, stimuli_info.positions', dotsize, black, [], []);
    drawCross(window1,W,H,barColor);
    tFixation = Screen('Flip', window1,[],1); % use 1 here to allow addition of the location!
    Screen('Flip', window1, tFixation + (fixationDuration + rand(1))- slack,1); %Adding between 0 and 1 second randomly to the minimum for the fixation duration
    
    % ----------- Prepare CUE
    drawCross(window1,W,H,barColor);
    Screen('FillOval', window1, barColorCue,[ (stimuli_info.positions(sequence_retrieval(t,1),1)-2*dotsize) (stimuli_info.positions(sequence_retrieval(t,1),2)-2*dotsize) (stimuli_info.positions(sequence_retrieval(t,1),1)+2*dotsize) (stimuli_info.positions(sequence_retrieval(t,1),2)+2*dotsize)]);
    
    % ----------- Send Trigger
    TriggerStart = GetSecs - ExpStart;
    sendTrigger(trg_type, trg_handle);
    
    % ----------- Show CUE
    % start of cue period
    tdum = Screen('Flip', window1); % Start of trial
    CueStart = tdum - ExpStart;
    if switch_behavior && sequence_retrieval(t,2) > 0
        CatchStart = CueStart;
    end
    
    if switch_standard
        % ----------- Reinstatement response
        waitForResponse = true;
        while waitForResponse %|| GetSecs - (CueStart+ExpStart) < cueTimeout+rand(1)
            [~,tdum,keyCode] = KbCheck;
            RetReinstStart = tdum - ExpStart;
            
            % Check for response keys
            if keyCode(KbName(responseKeys{1})) % remembered
                ret_reinst_resp = 1;
                RT_ret_reinst = RetReinstStart - CueStart;
                waitForResponse = false;
            elseif keyCode(KbName(responseKeys{2}))% forgotten
                ret_reinst_resp = 0;
                RT_ret_reinst = RetReinstStart - CueStart;
                waitForResponse = false;
            elseif keyCode(KbName(stopKey{1}))
                quitTask; return
            end
        end
        
        % Fill up a bar while subject keeps the image in mind
        centeredRect = CenterRectOnPointd([0 0 barWidth barHeight], xCenter,catch_y);
        growingsize = barWidth/(hz*retrievalTimeout);
        
        for fill = 1:hz*retrievalTimeout % fill up the bar
            Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
            drawCross(window1,W,H,barColor);
            Screen('FillRect', window1, [225 225 225], centeredRect);
            centeredRectFill = CenterRectOnPointd([0 0 growingsize*fill barHeight], xCenter,catch_y);
            Screen('FillRect', window1, [200 200 200], centeredRectFill);
            %         Screen('FillRect', window1, [200 200 200], [xCenter-barWidth/2 (H/1.3)-barHeight/2 xCenter-barWidth/2+growingsize*fill (H/1.3)+barHeight/2]);%centeredRectfill);
            Screen('Flip', window1,1)
        end
    elseif switch_behavior
        RT_ret_reinst = RetReinstStart - CueStart;
    end
    
    if sequence_retrieval(t,2) > 0 && switch_standard
        % ask catch question after reinstatement
        
        % show labels
        Screen('TextSize', window1, fontsize);
        if DPIswitch
            Screen('DrawText', window1, labels{1},(rand_positions(1,1)) - xoffset(1)+labeloffset, catch_y, textColor);%CORRECT
            Screen('DrawText',window1, '(forgotten)', W/2-100, catch_y+fontsize, [0 0 0]);
            Screen('DrawText', window1, labels{2},(rand_positions(2,1)) - xoffset(2)+labeloffset, catch_y, textColor);%LURE
        else
            DrawFormattedText_mod(window1, labels{1},'center', catch_y,textColor,[],[],[],[],[],[],(rand_positions(1,1)));
            DrawFormattedText_mod(window1, labels{2},'center', catch_y,textColor,[],[],[],[],[],[],(rand_positions(2,1)));
            DrawFormattedText(window1, '(forgotten)','center', catch_y+2*fontsize,textColor,[],[],[],[],[],[]);
        end
        
        tdum = Screen('Flip', window1);
        CatchStart = tdum - ExpStart;
    end
    
    if sequence_retrieval(t,2) > 0
        % wait for response
        waitForResponse = true;
        while waitForResponse
            [keyIsDown,respTime,keyCode] = KbCheck;
            pressedKeys = find(keyCode);
            
            if ~isempty(pressedKeys) && sum(ismember(KbName(responseKeys(2:end)),pressedKeys(1)))
                catch_resp_key = responseKeys{find(ismember(KbName(responseKeys(2:end)),pressedKeys(1)))+1};
                CatchRespStart = respTime - ExpStart;
                RT_catch = CatchRespStart - CatchStart;
                waitForResponse = false;
            elseif keyCode(KbName(stopKey{1}))
                quitTask; return
            end
        end
        
        % check answer
        catch_resp = 2; % no answer
        if rand_positions(1,1)==catch_positions(1,1) && catch_resp_key(1)==responseKeys{3}(1) %'L'
            %             correct_responses_catch = correct_responses_catch + 1;
            catch_resp = 1; %correct
        elseif rand_positions(1,1)==catch_positions(1,1) && catch_resp_key(1)==responseKeys{4}(1) %'R'
            catch_resp = 0; %incorrect
        elseif rand_positions(1,1)==catch_positions(2,1) && catch_resp_key(1)==responseKeys{4}(1) %'R'
            %             correct_responses_catch = correct_responses_catch + 1;
            catch_resp = 1;
        elseif rand_positions(1,1)==catch_positions(2,1) && catch_resp_key(1)==responseKeys{3}(1) %'L'
            catch_resp = 0;
        elseif catch_resp_key(1)==responseKeys{2}(1) %'D'
            catch_resp = 3; % forgotten
        end
        
        % --------- inter trial interval
        Screen('FillOval', window1, [255,255,255], [W/2-R H/2-R W/2+R H/2+R], []);
        ITIStart = Screen('Flip', window1);
        
        while GetSecs - ITIStart < postobjectTimeout
            [keyIsDown,secs,keyCode] = KbCheck;
            
            if keyCode(KbName(stopKey{1})) % quit the experiment
                quitTask; return
            end
        end
    end
    
    % ---------- print everything to log file
    stimID = sequence_retrieval(t,1);
    
    fprintf(outputfile, ['%s\t %d\t %s\t %d\t %d\t %d\t %f\t %f\t',...
        '%s\t %d\t %s\t %d\t %d\t ',...
        '%d\t %d\t %d\t %d\t %f\t %d\t %f\t %d\t %f\t',...
        '%d\t %f\t %s\t %f\t %d\t',...
        '%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t',...
        '%f\t %f\t %f\t %f\t %f\t %f\n'] ,...
        subID, sessionID, block_state, t, sequence_retrieval(t,2), stim_cue_idx(stimID), stimuli_info.positions(stimID,1), stimuli_info.positions(stimID,2),...
        stimuli_info.item_name{stimID}, stimuli_info.item_ID(stimID,:), imgList{stimID}, sum(sequence_retrieval(1:t,1)==stimID), stimuli_info.cat_exemplar(stimID),...
        stimuli_info.cat_perc_1(stimID), stimuli_info.cat_perc_2(stimID), stimuli_info.cat_sem_1(stimID),stimuli_info.cat_sem_2(stimID),...
        [], [], [], [], [],...
        ret_reinst_resp,RT_ret_reinst,catch_type,RT_catch,catch_resp,...
        SessionStart, [],[],RetrievalStart,TrialStart,CueStart,TriggerStart,[],...
        [],[],[], RetReinstStart,CatchStart,CatchRespStart);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% End the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DPIswitch
    Screen('DrawText',window1, 'Well done!',W/2-2*50,H/4, textColor);
else
    DrawFormattedText(window1, 'Well done!','center',H/4, textColor);
end
Screen('Flip',window1,1,[]);

KbStrokeWait;
KbReleaseWait();
RestrictKeysForKbCheck([]);
Screen(window1,'Close');
fclose('all');
close all
sca;
return

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Draw a fixation cross (overlapping horizontal and vertical bar)
function drawCross(window,W,H,barColor)

if nargin==3
    barColor = [0 0 0];
end

barLength = 16; % in pixels
barWidth = 2; % in pixels

Screen('FillRect', window, barColor,[ (W-barLength)/2 (H-barWidth)/2 (W+barLength)/2 (H+barWidth)/2]);
Screen('FillRect', window, barColor,[ (W-barWidth)/2 (H-barLength)/2 (W+barWidth)/2 (H+barLength)/2]);

end


function quitTask

KbReleaseWait();
fclose('all');
close all
sca

end