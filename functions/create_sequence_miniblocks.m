function [grand_sequence] = create_sequence_miniblocks...
    (nstim, nrep, mindiff, ntrials_block, nrep_Strials,mindiffS, skipblock)

% nstim = number of stimuli
% nrep = repetitions per stimuli
% prob = probability of category change

nblocks = (nstim*nrep) / ntrials_block;

mindist = min(nstim-1,mindiff); % number of trials the same stimulus should be apart

if mod(nstim*nrep, ntrials_block) ~= 0
    error('Cannot divide the trials into miniblocks')
end

%% make sequences per mini block

grand_sequence = zeros(nstim*nrep,2);

for b = 1:nblocks
    
    sequence = zeros(ntrials_block,2);
    
%     while length(unique(sequence(size(sequence,1)-mindist:size(sequence,1),1)))<=mindist %controlling that the same items is not repeated in the last 5 trials
        
        % create list of stimulus elements
        elements(:,1)               = 1:nstim;
        elements(:,2)               = 0; %counter;
        
        % all stimuli that need to be presented in this block
        all_elements = repmat(elements,ceil(ntrials_block/nstim),1);
        
        for i=1:size(all_elements,1)
            
            if i==1 % pick presentation slot for first element at random
                selected=randi(ntrials_block);
                sequence(i,:)=all_elements(selected,:);
                all_elements(selected,2)=1;
            else
                
                %Controlling that the next item is not repeated in n-1 and n-2 and
                %it belong to the category of the current trial
                if i<mindist+1
                    found=find(all_elements(:,2)~=1 & all_elements(:,1)~=sequence(i-1,1));
                else
                    dum = all_elements(:,2)~=1;
                    for d = 1:mindist
                       dum(all_elements(:,1)==sequence(i-d,1)) = 0;
                    end
                    found = find(dum);
%                     found=find(all_elements(:,2)~=1 & all_elements(:,1)~=sequence(i-1,1) & all_elements(:,1)~=sequence(i-2,1));
                end
                
                if ~isempty(found)
                    selected=found(randi(size(found,1)));
                    sequence(i,:)=all_elements(selected,:);
                    all_elements(selected,2)=1;
                else
                    found=find(all_elements(:,2)~=1);
                    selected=found(randi(size(found,1)));
                    sequence(i,:)=all_elements(selected,:);
                    all_elements(selected,2)=1;
                end
                
            end 
        end
%     end
    
    grand_sequence((b-1)*ntrials_block+1:b*ntrials_block,:) = sequence(1:ntrials_block,:);
end

%% determine which trials are special trials
% this section iteratively chooses the 'special' trials for each stimulus
% if it doesn't meet the criteria of spacing, it will try again, without
% warning.

if skipblock && nblocks==1
    error('Cannot asign drag-and-drop trials: not enough trials to skip a block')
end

contswitch = true;

while contswitch

    grand_sequence(:,2) = 0;
    
    for s = 1:nstim
        
        % find stimuli
        if skipblock
            id = find(grand_sequence(ntrials_block+1:end,1) == s) + ntrials_block;
        else
            id = find(grand_sequence(1:end,1) == s);
        end
        goodid = [];
        
        % test all possible trials for this stimulus for spacing
        incl = [];
        for t = 1:length(id)
            tid = id(t);
            incl(t) = sum(grand_sequence(tid-mindiffS:min(length(grand_sequence),tid+mindiffS),2)) == 0;
        end
        
        % if there are not enough trials to meet spacing, break and try
        % again
        if sum(incl) < nrep_Strials
            break
        end
        
        % find remaining ids and choose at random
        dum = id(find(incl));
        goodid = dum(randperm(length(dum),nrep_Strials));
       
        % store
        grand_sequence(goodid,2) = 1;
    end
    
    if s == nstim && length(goodid) == nrep_Strials
        contswitch = false;
    end
end

end