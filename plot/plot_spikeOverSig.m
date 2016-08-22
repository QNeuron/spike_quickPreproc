function plot_spikeOverSig(spikeMat,Channel,thresholds,TimeLimit,expInfo,dataPath)

sr = expInfo.exp.dataDeviceSampleRate;
binDur = 1/sr;
abscissa = TimeLimit(1):binDur:TimeLimit(2);

files = excludeDots(dir(dataPath));
datFiles = {files(arrayfun(@(x)(~isempty(strfind(x.name,'.f32.dat'))),files)).name};

filterLFP = designfilt('lowpassiir', ...        % Response type
    'PassbandFrequency',500, ...     % Frequency constraints
    'StopbandFrequency',500+150, ...
    'PassbandRipple',4, ...          % Magnitude constraints
    'StopbandAttenuation',65, ...
    'DesignMethod','butter', ...      % Design method
    'MatchExactly','passband', ...   % Design method options
    'SampleRate',sr);               % Sample rate

filterAP = designfilt('bandpassiir', ...       % Response type
    'StopbandFrequency1',500, ...    % Frequency constraints
    'PassbandFrequency1',500+150, ...
    'PassbandFrequency2',10000, ...
    'StopbandFrequency2',11000, ...
    'StopbandAttenuation1',65, ...   % Magnitude constraints
    'PassbandRipple',1, ...
    'StopbandAttenuation2',65, ...
    'DesignMethod','butter', ...      % Design method
    'MatchExactly','passband', ...   % Design method options
    'SampleRate',sr);               % Sample rate

sweepDur = expInfo.sweepLength/sr;
NSweep = length(datFiles);
lengthVect = (1:NSweep)*sweepDur;
SweepToLoad = find(lengthVect>=TimeLimit(1)&lengthVect<=TimeLimit(2)+sweepDur);

zed = regexp(datFiles{1},'\.','split');
SweepInd = find(strcmp(zed,'sweep'))+1;
for i = 1:NSweep,
    zed = regexp(datFiles{i},'\.','split');
    SweepNb(i) = str2double(zed{SweepInd});
end

s = [];
for i = SweepToLoad,
    dataS = f32read(fullfile(dataPath,datFiles{SweepNb==i}));
    dataS = dataS(Channel:32:end);
    s = [s; dataS];
end

sAP = filter(filterAP,s);
sLFP = filter(filterLFP,s);

realAbsc = (lengthVect(SweepToLoad(1))-sweepDur):binDur:(lengthVect(SweepToLoad(end))-binDur);

spkTimes = spikeMat(spikeMat(:,5)==Channel & spikeMat(:,1)>=abscissa(1) & spikeMat(:,1)<=abscissa(end),1);

figure
plot(abscissa,sAP(realAbsc>=abscissa(1)&realAbsc<=abscissa(end)+binDur));
hold on
plot(spkTimes,ones(1,length(spkTimes))*max(sAP),'color','r','linestyle','none','marker','*');
line([abscissa(1) abscissa(end)],[thresholds(Channel) thresholds(Channel)],'color','k')
hold off
title('MUA signal')


figure
plot(abscissa,sLFP(realAbsc>=abscissa(1)&realAbsc<=abscissa(end)+binDur));
hold on
plot(spkTimes,ones(1,length(spkTimes))*max(sLFP),'color','r','linestyle','none','marker','*');
hold off
title('LFP signal')


figure
plot(abscissa,s(realAbsc>=abscissa(1)&realAbsc<=abscissa(end)+binDur));
hold on
plot(spkTimes,ones(1,length(spkTimes))*max(s),'color','r','linestyle','none','marker','*');
hold off
title('Raw signal')


end