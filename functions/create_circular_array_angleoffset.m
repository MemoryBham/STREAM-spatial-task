
function [cue_positions_sort,cue_positions_rand] = create_circular_array_angleoffset(x,y,r,a_offset,stimuli)

nitems = size(stimuli,1);

% Making the circle coordinates
th = 0:pi/200:2*pi;
xunit = r * cos(th+a_offset) + x;
yunit = r * sin(th+a_offset) + y;
xyunits=[xunit; yunit]';

% Create an equal distance between all possible items' positions
distance=round(400/nitems);
ao=1:distance:400;
cue_positions_sort=xyunits(ao,:); %List of position for the items

% randomly assign positions
stim_cue_idx = randperm(size(cue_positions_sort,1));
cue_positions_rand = cue_positions_sort(stim_cue_idx,:);

if sum(size(stimuli)>1)>1
% make sure left-facing items are not all on the left etc.
leftpos = find(cue_positions_rand(:,1)<x); % left half of the circle
rightpos = find(cue_positions_rand(:,1)>x); % right half of the circle

if sum(stimuli.cat_perc_2(leftpos)==1) == length(leftpos)  ||  ...% left-facing stimuli at left positions
        sum(stimuli.cat_perc_2(rightpos)==2) == length(rightpos) % right-facing stimuli at right positions
    % swap two locations at random
    lid = find(stimuli.cat_perc_2(leftpos)==1); 
    rid = find(stimuli.cat_perc_2(rightpos)==2);
    if ~isempty(lid) && ~isempty(rid)
        idl = leftpos(lid(randperm(length(lid),1)));
        idr = rightpos(rid(randperm(length(rid),1)));
        dum = cue_positions_rand(idr,:);
        cue_positions_rand(idr,:) = cue_positions_rand(idl,:);
        cue_positions_rand(idl,:) = dum;
    end
end
end



 