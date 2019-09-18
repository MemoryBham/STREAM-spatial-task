function stimuli = create_stimulus_split(stimuli,nsessions, nstim)

% for 4 stimuli:
% 1. random first choice
% 2. exclude for this session:
    % - same combination of perc1 & perc2
    % - same sem1/sem2 with same perc1 & perc2 combination
% 3. exclude for other session:
    % - same cat. with perc1-perc2 combination in other session
    
% there are 4 possible combinations of perc1 & perc2    

% for 8 stimuli:
% 1. random first choice
% 2. exclude for this session:
    % - same perc1/perc2 for this exemplar category within this session
    % - same sem1/sem2 with same perc1 & perc2 combination
% 3. exclude for other session:
    % - same cat. with perc1-perc2 combination in other session

% for 16 stimuli: 
% 1. random first choice
% 2. exclude for this session:
    % - same combination of perc1 & perc2 for this category
% 3. exclude for other session:
    % - same cat. with perc1-perc2 combination in other session
    
stimuli.session = zeros(size(stimuli,1),1);
available_tot = ones(size(stimuli,1),1);

% loop over sessions
for s = 1:nsessions
    
    % available stimuli for this category
    available_ses = stimuli.session==0 & available_tot;
    
    % within exemplar category
    for ce = unique(stimuli.cat_exemplar)'
        
        % exclude other exemplar categories
        available_cat = available_ses;
        available_cat(stimuli.cat_exemplar~=ce) = 0;
%         sum(available_cat)
        
        for st = 1:nstim/4
        
            % select stimulus at random from remaining available options
            adum = find(available_cat);
            tid = adum(randperm(length(adum),1));
            stimuli.session(tid) = s;
            
            % exclude WITHIN EXEMPLAR category, WITHIN session
            % exclude same exemplar from session
            id_e = strcmp(stimuli.label_exemplar,stimuli.label_exemplar(tid)); 
            available_cat(id_e) = 0;
            % exclude same perc1&perc2 combination within exemplar
            % category, within this session
            id_epp = stimuli.cat_exemplar==stimuli.cat_exemplar(tid) & ...
                stimuli.cat_perc_1==stimuli.cat_perc_1(tid) & stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            available_cat(id_epp) = 0;
            
            % exlude same perc1 or perc2 within exemplar category, within
            % session if a threshold is reached
            id_ep1 = stimuli.cat_exemplar==stimuli.cat_exemplar(tid) & stimuli.cat_perc_1==stimuli.cat_perc_1(tid); 
            if sum(stimuli.session(id_ep1)==s) >= nstim/8
                available_cat(id_ep1) = 0;
            end   
            id_ep2 = stimuli.cat_exemplar==stimuli.cat_exemplar(tid) & stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            if sum(stimuli.session(id_ep2)==s) >= nstim/8
                available_cat(id_ep2) = 0;
            end
            
            % exclude ACROSS EXEMPLAR categories, WITHIN session
            % exclude same combination of perc1 & perc2 for this session for 
            % all exemplar categories if a threshold is reached (4 or 8 stimuli) 
            id_pp = stimuli.cat_perc_1==stimuli.cat_perc_1(tid) & stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            if sum(stimuli.session(id_pp)==s) >= nstim/4
                % exclude all other pairs of perc1 & perc2 for this session
                available_ses(id_pp) = 0;
            end
            
            % exclude same perc1/perc2 from the same session if a threshold 
            % is reached
            id_p1 = stimuli.cat_perc_1==stimuli.cat_perc_1(tid); 
            if sum(stimuli.session(id_p1)==s) >= nstim/2
                available_ses(id_p1) = 0;
            end
            id_p2 = stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            if sum(stimuli.session(id_p2)==s) >= nstim/2
                available_ses(id_p2) = 0;
            end
            
            % exclude WITHIN SEMANTIC categories, WITHIN session
            if nstim>4
            % exclude same combination of sem1/sem2 & perc1/perc2 within
            % session if a threshold is reached
            id_s1p1 = stimuli.cat_sem_1==stimuli.cat_sem_1(tid) & stimuli.cat_perc_1==stimuli.cat_perc_1(tid); 
            if sum(stimuli.session(id_s1p1)==s) >= nstim/4
                available_ses(id_s1p1) = 0;
            end
            id_s1p2 = stimuli.cat_sem_1==stimuli.cat_sem_1(tid) & stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            if sum(stimuli.session(id_s1p2)==s) >= nstim/4
                available_ses(id_s1p2) = 0;
            end
            id_s2p1 = stimuli.cat_sem_2==stimuli.cat_sem_2(tid) & stimuli.cat_perc_1==stimuli.cat_perc_1(tid); 
            if sum(stimuli.session(id_s2p1)==s) >= nstim/4
                available_ses(id_s2p1) = 0;
            end
            id_s2p2 = stimuli.cat_sem_2==stimuli.cat_sem_2(tid) & stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            if sum(stimuli.session(id_s2p2)==s) >= nstim/4
                available_ses(id_s2p2) = 0;
            end
            end
            
            id_s1p1p2 = stimuli.cat_sem_1==stimuli.cat_sem_1(tid) & ...
                stimuli.cat_perc_1==stimuli.cat_perc_1(tid) & stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            if sum(stimuli.session(id_s1p1p2)==s) >= nstim/8
                available_ses(id_s1p1p2) = 0;
            end
            id_s2p1p2 = stimuli.cat_sem_2==stimuli.cat_sem_2(tid) & ...
                stimuli.cat_perc_1==stimuli.cat_perc_1(tid) & stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            if sum(stimuli.session(id_s2p1p2)==s) >= nstim/8
                available_ses(id_s2p1p2) = 0;
            end
            
            % exclude WIHTIN EXEMPLAR category, ACROSS sessions
            % exclude same exemplar from all sessions (only 1 session
            % possible for 16 stimuli)
            id_e = strcmp(stimuli.label_exemplar,stimuli.label_exemplar(tid)); 
            available_tot(id_e) = 0;
            % if needed: exclude same combination of perc1 & perc2 within exemplar 
            % category for other sessions (4 or 8 stimuli)
            id_epp = stimuli.cat_exemplar==stimuli.cat_exemplar(tid) & ...
                stimuli.cat_perc_1==stimuli.cat_perc_1(tid) & stimuli.cat_perc_2==stimuli.cat_perc_2(tid); 
            if sum(stimuli.session(id_epp)>0) >= nstim/8
                available_tot(id_epp) = 0;
            end
        end
        
    end
end
    
    
% %% old version
% % CAUTION: this code only works for 4 or 8 stimuli!
% 
% % for two sessions: 
% % CRITERIA
% % two of each exemplar cat
% % one left, one right per exemplar cat
% % one photo, one drawing per exemplar cat
% 
% stimuli.session = zeros(size(stimuli,1),1);
% 
% % loop over sessions
% for s = 1:nsessions
%     
%     % within exemplar category
%     for ce = unique(stimuli.cat_exemplar)'
%         
%         % available stimuli for this category
%         available = stimuli.cat_exemplar==ce & stimuli.session==0;
%                
%         % check balancing of perceptual categories across exemplar
%         % categories
%         if nstim == size(stimuli,1)/2 && sum(stimuli.cat_perc_1(stimuli.session==s)==1)==nstim/2
%             available(stimuli.cat_perc_1==1) = 0;
%         elseif nstim == size(stimuli,1)/2 && sum(stimuli.cat_perc_1(stimuli.session==s)==2)==nstim/2
%             available(stimuli.cat_perc_1==2) = 0;
%         end
%         % check balancing of combinations of perceptual categories within
%         % session
%         if sum(stimuli.cat_perc_1(stimuli.session==s)==1 & stimuli.cat_perc_2(stimuli.session==s)==1)==nstim/4
%             available(stimuli.cat_perc_1==1 & stimuli.cat_perc_2==1) = 0;
%         elseif sum(stimuli.cat_perc_1(stimuli.session==s)==1 & stimuli.cat_perc_2(stimuli.session==s)==2)==nstim/4
%             available(stimuli.cat_perc_1==1 & stimuli.cat_perc_2==2) = 0;
%         end
%         if sum(stimuli.cat_perc_1(stimuli.session==s)==2 & stimuli.cat_perc_2(stimuli.session==s)==1)==nstim/4
%             available(stimuli.cat_perc_1==2 & stimuli.cat_perc_2==1) = 0;
%         elseif sum(stimuli.cat_perc_1(stimuli.session==s)==2 & stimuli.cat_perc_2(stimuli.session==s)==2)==nstim/4
%             available(stimuli.cat_perc_1==2 & stimuli.cat_perc_2==2) = 0;
%         end
%         
%         % check balancing of combinations of perceptual categories across
%         % session
%         if sum(stimuli.cat_perc_1(available==1)==1 & stimuli.cat_perc_2(available==1)==1)==nstim/4-1
%             available(stimuli.cat_perc_1==1 & stimuli.cat_perc_2==1) = 0;
%         elseif sum(stimuli.cat_perc_1(available==1)==1 & stimuli.cat_perc_2(available==1)==2)==nstim/4-1
%             available(stimuli.cat_perc_1==1 & stimuli.cat_perc_2==2) = 0;
%         end
%         if sum(stimuli.cat_perc_1(available==1)==2 & stimuli.cat_perc_2(available==1)==1)==nstim/4-1
%             available(stimuli.cat_perc_1==2 & stimuli.cat_perc_2==1) = 0;
%         elseif sum(stimuli.cat_perc_1(available==1)==2 & stimuli.cat_perc_2(available==1)==2)==nstim/4-1
%             available(stimuli.cat_perc_1==2 & stimuli.cat_perc_2==2) = 0;
%         end
%         
%         % start with a random stimulus
%         adum = find(available);
%         tid = adum(randperm(length(adum),1));
%         stimuli.session(tid) = s;
%         
%         available(strcmp(stimuli.label_exemplar,stimuli.label_exemplar(tid))) = 0;
%         
%         % for >4 stimuli
%         if nstim > length(unique(stimuli.cat_exemplar))
%             % update stimulus availability
%             if stimuli.cat_perc_1(tid) == 1
%                 available(stimuli.cat_perc_1==1) = 0;
%             else
%                 available(stimuli.cat_perc_1==2) = 0;
%             end
%             if stimuli.cat_perc_2(tid) == 1
%                 available(stimuli.cat_perc_2==1) = 0;
%             else
%                 available(stimuli.cat_perc_2==2) = 0;
%             end
%             
%             % find second stimulus
%             adum = find(available);
%             tid = adum(randperm(length(adum),1));
%             stimuli.session(tid) = s;
%         end
%         
%     end
% end


end