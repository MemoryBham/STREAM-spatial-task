
imageFolder = './stimuli/stimuli_dogs_birds_cars_planes/';

nsessions = 2;
nstim = 8;
nsub = 12;
nrep = 12;

stim_store = [];

for i = 1:nrep
    
    stimuli_info = readtable(strcat(imageFolder, '/stimuli_info.txt'));
    % create random stim split
    stimuli_info = create_stimulus_split(stimuli_info,nsessions, nstim);
    
    % store 
    stim_store(:,i) = stimuli_info.session;
end

figure; imagesc(stim_store)

% per stimulus
figure; plot(sum(stim_store>0,2)); ylim([0,20])

% per exemplar
c = 1;
for e = unique(stimuli_info.label_exemplar)'
   perc1_1(c) = sum( sum(stim_store>0,2) .* double(stimuli_info.cat_perc_1==1) .* strcmp(stimuli_info.label_exemplar,e));
   perc1_2(c) = sum( sum(stim_store>0,2) .* double(stimuli_info.cat_perc_1==2) .* strcmp(stimuli_info.label_exemplar,e));
   perc2_1(c) = sum( sum(stim_store>0,2) .* double(stimuli_info.cat_perc_2==1) .* strcmp(stimuli_info.label_exemplar,e));
   perc2_2(c) = sum( sum(stim_store>0,2) .* double(stimuli_info.cat_perc_2==2) .* strcmp(stimuli_info.label_exemplar,e));
   c = c+1;
end

figure; 
subplot(211); hold on
plot(perc1_1)
plot(perc1_2)
title('Drawing - photograph')
subplot(212); hold on
plot(perc2_1)
plot(perc2_2)
title('Left - right')

% per session
[pair12,pair13,pair14,pair23,pair24,pair34] = deal(zeros(4,1));
for c = 1:4
    ex = unique(stimuli_info.label_exemplar(stimuli_info.cat_exemplar==c));
    for s = 1:nsessions
       dum1 = double(stim_store==s) .* double(repmat(strcmp(stimuli_info.label_exemplar,ex(1)),[1,nrep]));
       dum2 = double(stim_store==s) .* double(repmat(strcmp(stimuli_info.label_exemplar,ex(2)),[1,nrep]));
       dum3 = double(stim_store==s) .* double(repmat(strcmp(stimuli_info.label_exemplar,ex(3)),[1,nrep]));
       dum4 = double(stim_store==s) .* double(repmat(strcmp(stimuli_info.label_exemplar,ex(4)),[1,nrep]));
       
       pair12(c,:) = pair12(c,:)+sum(sum(dum1,1) .* sum(dum2,1));
       pair13(c,:) = pair13(c,:)+sum(sum(dum1,1) .* sum(dum3,1));
       pair14(c,:) = pair14(c,:)+sum(sum(dum1,1) .* sum(dum4,1));
       pair23(c,:) = pair23(c,:)+sum(sum(dum2,1) .* sum(dum3,1));
       pair24(c,:) = pair24(c,:)+sum(sum(dum2,1) .* sum(dum4,1));
       pair34(c,:) = pair34(c,:)+sum(sum(dum3,1) .* sum(dum4,1));
    end
end

figure;
bar([pair12,pair13,pair14,pair23,pair24,pair34])