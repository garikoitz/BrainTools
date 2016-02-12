%% Lexicality experiment
%
% Here we are manipulating two stimulus properties (a) image contrast and
% (b) frequency of the letter combinations making up the text. Contrast
% effects the perceptual decision while more frequent letter combinations
% effect the lexical decision. Stimuli will be presented with two tasks:
% (1) fixation task and (2) lexical decision task (word or pseudoword). We
% will have 3 levels of contrast (low, medium, high) and 4 levels of letter
% frequency (consonant strings, infrequent pseudowords, frequent
% pseudowords, high frequency real words).
% This makes for 12 total stimulus categories (3 contrast X 4 lexicality).
% Each stimulus category will be shown 4 times per run.
% and two tasks.

nstim = 120 %The total number of stimuli needed for each category
sorttype = 'UN3_F';% The header category to sort based on

% Load in excel spreadsheet (high freq words)
[n,t,c]=xlsread('Lexicality_Stimuli.xlsx',1);
% Header
h = t(1,:);
% letter string
lstr = t(2:end,1);
% Sort the letter strings based on their orthographic frequency
[~, oi] = sort(n(:,find(strcmp(sorttype,h))-1));
% Put the stimuli into a cell array
st(:,1) = lstr(oi(end-nstim+1:end));

%% Load in excel spreadsheet (trigrams)
[n,t,c]=xlsread('Lexicality_Stimuli.xlsx',3);
% Header
h = t(1,:);
% letter string
lstr = t(2:end,1);
% Sort the letter strings based on their orthographic frequency
[~, oi] = sort(n(:,find(strcmp(sorttype,h))-1));
% Put the high frequency trigrams into a cell
st(:,2) = lstr(oi(end-nstim+1:end));
% And the low frequency trigrams into another cell
st(:,3) = lstr(oi(1:nstim));

%% Load in excel spreadsheet (bigrams)
[n,t,c]=xlsread('Lexicality_Stimuli.xlsx',5);
% Header
h = t(1,:);
% letter string
lstr = t(2:end,1);
% Sort the letter strings based on their orthographic frequency
[~, oi] = sort(n(:,find(strcmp(sorttype,h))-1));
% Add the low frequency bigrams into another cell
st(:,3) = lstr(oi(1:nstim));

%% Load in excel spreadsheet (cons str)
[n,t,c]=xlsread('Lexicality_Stimuli.xlsx',6);
lstr = t(2:end,1);
% grab random cons str
st(:,4) = lstr(randsample(1:length(lstr),nstim));

%% Now shuffle the rows
for c = 1:4
    r = Shuffle(1:size(st,1));
    st(:,c) = st(r,c);
end
%% render the stimuli
clevels = [255 136 131]
bg = 128;
stimcon(1:40, 1:4)   = 1;
stimcon(41:80, 1:4)  = 2;
stimcon(81:120, 1:4) = 3;
c = 0; % counter
for ii = 1:size(st,1)
    for jj = 1:size(st,2)
        c = c+1;
        tmp = uint8(renderText(st{ii,jj},'Courier',24,4));
        % Set contrast
        tmp(tmp==1)=clevels(stimcon(ii,jj));
        tmp(tmp==0)=bg;
        % stack it up
        img(:,:,c) = tmp;
        sc(c) = stimcon(ii,jj);
        sl(c) = jj;
        stimcat(c) = (jj-1).*4+stimcon(ii,jj);
    end
end
stimcat = sl + (sc-1).*4

%% Decide on the frame order
nreps = 5; %number of times that each stimulus category will be repeated in a run
nruns = 8; %number of runs to make
nblank = 10; %number of blanks in a run
cnums = unique(stimcat); %all the stimulus category numbers
for rn = 1:nruns
    tmp = [];
    for c = cnums
       cidx = find(stimcat == c)
       tmp = horzcat(tmp, cidx(rn.*nreps-4:rn.*nreps));
    end
    % Add in blanks and Shuffle stim order to randomize conditions for the run
    stimorder(rn,:) = Shuffle(horzcat(tmp, zeros(1,nblank)));
end

%% Save it out
desc = ' img is the image stack \n sc denotes the contrast of each image\n sl denotes the lexicality level\n stimcat denotes the category\n stimorder gives the order of the stimulus for each run\n'
save LexicalityExp img sc sl stimcat desc stimorder

