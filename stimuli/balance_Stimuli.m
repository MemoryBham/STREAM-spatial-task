% to balance:

% is stimulus left/right or photograph/drawing
% balanced every 4 subjects

% is exemplar paired with what other examplar?
% balanced every 3 subjects (6 to correct for session ID)

% total stimulus set is balanced within 16 subjects

% is stimulus located in location 1-16?
% balanced every 16 subjects

% this code is for:
% # stimuli 16
% spread over 2 sessions

imageFolder = './stimuli/stimuli_dogs_birds_cars_planes/';
stimuli_info = readtable(strcat(imageFolder, '/stimuli_info.txt'));

%% make session decision matrix - 2 sessions

% the base unit: all possible pairs of 2 exemplars
cb = nchoosek(4,2);
Dses_base = ones(4, cb);
dum = combnk(1:4,2);
for i = 1:size(dum,1)
    Dses_base(dum(i,:),i) = 2;
end
tmp = Dses_base(:,3);
Dses_base(:,3) = Dses_base(:,1);
Dses_base(:,1) = tmp;
Dses_base = repmat(Dses_base,[1,3]);

% for loops over the 4 catergories to make complete decision matrix
Dses = [];%zeros(16,cb^4);
for j = 1:cb^3
    
   cp2 = mod(j,cb);
   cp3 = mod(j+1,cb);
   cp4 = mod(j+2,cb);
    
   cc1 = mod(j,cb) + 1; 
   cc2 = mod(j,cb) + 1 + cp2;
   cc3 = mod(j,cb) + 1 + cp3;
   cc4 = mod(j,cb) + 1 + cp4;
   
   Dtemp = cat(1, Dses_base(:,cc1:cc1+cb-1), Dses_base(:,cc2:cc2+cb-1), ...
       Dses_base(:,cc3:cc3+cb-1),Dses_base(:,cc4:cc4+cb-1));
   
   Dses = cat(2,Dses,Dtemp);
end

%% make perc decision matrix - 2 sessions

Dsem = [12,12,12,12,11,11,11,11,22,22,22,22,21,21,21,21]';

pr1 = [11,22];
pr2 = [12,21];
Dperc_base = [pr1,pr2];
Dperc_base_r = [fliplr(pr1),fliplr(pr2)];
Dperc_op = [pr2,pr1];
Dperc_op_r = [fliplr(pr2), fliplr(pr1)];

Dperc = zeros(size(Dses));
% start with the first row
stid = 1;
Dperc(stid,:) = repmat(Dperc_base,[1,size(Dperc,2)/size(Dperc_base,2)]);

for j = 1:size(Dperc,2)
    stid = 1;
    id = mod(j,4);
    if id == 0; id = 4; end
    
    % ---- for the session defined by row 1
    % this defines the perc of the exemplar from the same category in the same
    % session
    Dperc(stid+find(Dses(2:4,j)==Dses(stid,j)),j) = Dperc_base_r(id);

    % it also defined the perc of the stimuli with the same sem_1 and sem_2
    % same sem_1
    s1id = 4+find(round(Dsem(5:end)/10) == round(Dsem(1)/10) & Dses(5:end,j)==Dses(stid,j));
    Dperc(s1id(1),j) = Dperc_op_r(id);
    Dperc(s1id(2),j) = Dperc_op(id);
    
    % same sem_2
    s2id = 4+find(mod(Dsem(5:end),10) == mod(Dsem(1),10) & Dses(5:end,j)==Dses(stid,j));
    Dperc(s2id(1),j) = Dperc_op(id);
    Dperc(s2id(2),j) = Dperc_op_r(id);
    
    % both different sem_1 and sem_2
    snid = 4+find(round(Dsem(5:end)/10) ~= round(Dsem(1)/10) & ...
        mod(Dsem(5:end),10) ~= mod(Dsem(1),10) & Dses(5:end,j)==Dses(stid,j));
    Dperc(snid(1),j) = Dperc_base_r(id);
    Dperc(snid(2),j) = Dperc_base(id);
    
    % ---- now do first of the opposite session
    % find first row with Dperc==0
    stid = 1+find(Dperc(stid+1:stid+3,j)==0, 1, 'first');
    Dperc(stid,j) = Dperc_op_r(id);
    
    % fill in the last gap for the first category
    Dperc(Dperc(1:4,j)==0,j) = Dperc_op(id);
    
    % fill up the same semantic categories
    s1id = 4+find(round(Dsem(5:end)/10) == round(Dsem(1)/10) & Dses(5:end,j)==Dses(stid,j));
    Dperc(s1id(1),j) = Dperc_base(id);
    Dperc(s1id(2),j) = Dperc_base_r(id);
    
    % same sem_2
    s2id = 4+find(mod(Dsem(5:end),10) == mod(Dsem(1),10) & Dses(5:end,j)==Dses(stid,j));
    Dperc(s2id(1),j) = Dperc_base_r(id);
    Dperc(s2id(2),j) = Dperc_base(id);
    
    % both different sem_1 and sem_2
    snid = 4+find(round(Dsem(5:end)/10) ~= round(Dsem(1)/10) & ...
        mod(Dsem(5:end),10) ~= mod(Dsem(1),10) & Dses(5:end,j)==Dses(stid,j));
    Dperc(snid(1),j) = Dperc_op(id);
    Dperc(snid(2),j) = Dperc_op_r(id);
end


% for i = 1:4
%     stid = (i-1)*4+1;
%     Dperc(stid,:) = repmat(Dperc_base,[1,size(Dperc,2)/size(Dperc_base,2)]);
%     
%     for j = 1:size(Dperc,2)
%        id = mod(j,4);
%        if id == 0; id = 4; end
%         
%        % find Dperc(stid:stid+3,:) with same Dses as stid
%        Dperc(stid+find(Dses(stid+1:stid+3,j)==Dses(stid,j)),j) = Dperc_base_r(id);
%               
%        % find first row with Dperc==0
%        Dperc(stid+find(Dperc(stid+1:stid+3,j)==0, 1, 'first'),j) = Dperc_op(id);
%        
%        % fill up remaining spaces with 
%        Dperc(stid+find(Dperc(stid+1:stid+3,j)==0),j) = Dperc_op_r(id);
%     end
% end

%% make location decision matrix - 2 sessions
locs1 = 1:8;
locs2 = 9:16;
alllocs1 = perms(locs1)'; 
alllocs2 = perms(locs2)'; 

% to prevent right-facing always right and left-facing always left
rlocs = [1,2,8,9,10,15,16];
llocs = [4,5,6,11,12,13,14];

sdiff = 5; % max is 7, but this is unlikely to converge!

Dloc = zeros(size(Dses));

id1 = randi(size(alllocs1,2),1);
Dloc(Dses(:,1)==1,1) = alllocs1(:,id1);
alllocs1 = alllocs1(:,setdiff(1:size(alllocs1,2),id1));

id2 = randi(size(alllocs1,2),1);
Dloc(Dses(:,1)==2,1) = alllocs2(:,id2);
alllocs2 = alllocs2(:,setdiff(1:size(alllocs2,2),id1));


for i = 2:size(Dloc,2)
    
    for s = 1:2
        locs = eval(['locs',num2str(s)]);
        alllocs = eval(['alllocs', num2str(s)]);
        
        % find ids for this session
        ids = find(Dses(:,i) == s);
        
        % find distributions across locations for all stimuli
        Lhist = zeros(8);
        for l = 1:8
            Lhist(:,l) = sum(Dloc(ids,1:i-1)==locs(l),2);
        end
        
        % find the lowest and the highest
        [stim_min,loc_min] = ind2sub(size(Lhist),find(Lhist == min(Lhist(:)), 1));
        
        tmp = setdiff(1:8,stim_min);
        [~, tmp2] = max(sum(Lhist(tmp,:) == max(Lhist(:)),2));
        stim_max = tmp(tmp2);
        notloc_max = find(Lhist(stim_max,:) == max(Lhist(stim_max,:)));
        
        % find a location permutation that fills up the lowest and doesn't add
        % to the highest
        goodperms = find(alllocs(stim_min,:) == locs(loc_min));
        badperms = find(ismember(alllocs(stim_max,goodperms), locs(notloc_max)));
        goodperms = goodperms(setdiff(1:length(goodperms),badperms));
        
        % check that it meets criteria:
        
        % do not repeat locations of previous sdiff subjects
        if i > sdiff+1
        for t = 1:sdiff
            badperms = find(sum(alllocs(:,goodperms) == repmat(Dloc(ids,i-t), [1,length(goodperms)]),1)>0);
            goodperms = goodperms(setdiff(1:length(goodperms),badperms));
        end
        end
        
        % remove perms with all right-facing stim on the right
        % left- and right-facing stim
        lid = find(mod(Dperc(ids,i),10) == 1);
        rid = find(mod(Dperc(ids,i),10) == 2);
        badperms = find(sum(ismember(alllocs(lid,goodperms),llocs),1)>=3 | sum(ismember(alllocs(rid,goodperms),rlocs),1)>=3);
        goodperms = goodperms(setdiff(1:length(goodperms),badperms));
        
        % pick perm at random
        idl = goodperms(randi(length(goodperms),1));
        Dloc(ids,i) = alllocs(:,idl);
        
        % remove location permutation from list
        eval(['alllocs', num2str(s), ' = alllocs(:,setdiff(1:size(alllocs,2),idl));']);
        
    end
end

% test
% figure; hold on
% for i = 1:16
%     plot(cumsum(Dloc(1,:) == i))
% end
% plot(1:16:1296,1:(1296/16), 'k', 'linewidth',2)

%% convert everything to tables

stimuli_sessions = zeros(size(stimuli_info,1),size(Dloc,2));
stimuli_locations = zeros(size(stimuli_info,1),size(Dloc,2));
stimlabels = unique(stimuli_info.label_exemplar, 'stable');

for i = 1:size(Dloc,2)
   for j = 1:size(Dloc,1)
      % find stimID
      cat1 = round(Dperc(j,i)/10);
      cat2 = mod(Dperc(j,i),10);
      stimID = find(strcmp(stimuli_info.label_exemplar,stimlabels(j)) & ...
          stimuli_info.cat_perc_1 == cat1 & stimuli_info.cat_perc_2 == cat2);
      
      % store session for stimID
      stimuli_sessions(stimID,i) = Dses(j,i);
      
      % store location for stimID
      stimuli_locations(stimID,i) = Dloc(j,i);
   end
end

%% store locations and selection

save([imageFolder, '/stimuli_sessions.mat'], 'stimuli_sessions')
save([imageFolder, '/stimuli_locations.mat'], 'stimuli_locations')
