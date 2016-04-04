function data = get_spike_quick(dataPath)

%% Get data ID

if ~exist('dataPath','var') || isempty(dataPath),
    dataPath = cd;
end

files = excludeDots(dir(dataPath));
datFiles = {files(arrayfun(@(x)(~isempty(strfind(x.name,'.f32.dat'))),files)).name};


load(fullfile(dataPath,'gridInfo.mat'))
nChannels = expt.nChannels;
nSweep = grid.nSweepsDesired;
sr = expt.dataDeviceSampleRate;

%% LFP removal (low pass / high pass)
lfpFreqLim = 500; % Hz
plotDur = 10; %sec

filterLFP = designfilt('lowpassiir','FilterOrder',8, ...
    'PassbandFrequency',lfpFreqLim,'PassbandRipple',1, ...
    'SampleRate',sr);
filterAP = designfilt('bandpassfir','FilterOrder',20, ...
    'CutoffFrequency1',lfpFreqLim-100,'CutoffFrequency2',10000, ...
    'SampleRate',sr,'PassbandRipple',1);
% fvtool(filterLFP);

data = struct;

% Add condition on nSweep / nRep to avoid to long loading time.
dataS = f32read(fullfile(dataPath,datFiles{1})); % Load first sweep to get sweep length
for ch = 1:nChannels,
        dataCh = dataS(ch:nChannels:end); % de-interleave
        LFP = filter(filterLFP,dataCh);
    %         LFP = filtfilt(filterLFP,dataCh);
        AP = filter(filterAP,dataCh);
        
        data(ch).unfiltTrace{1} = dataCh;
        data(ch).LFP{1} = LFP;
        data(ch).AP{1} = AP;
end
sweepDur = length(dataS)/nChannels/sr; % seconds
nSweepNeeded = ceil(plotDur/sweepDur);

for i = 2:nSweepNeeded, % Load more seep if more are needed
    dataS = f32read(fullfile(dataPath,datFiles{i}));
    
    for ch = 1:nChannels,
        dataCh = dataS(ch:nChannels:end); % de-interleave
        LFP = filter(filterLFP,dataCh);
    %         LFP = filtfilt(filterLFP,dataCh);
        AP = filter(filterAP,dataCh);
        
        data(ch).unfiltTrace{i} = dataCh;
        data(ch).LFP{i} = LFP;
        data(ch).AP{i} = AP;
    end
    
    
end

defThr = 0;

handles.hF = figure('units','normalized','outerposition',[0 0 1 1]);
% [pos]=subplot_pos(nChannels,1,edgel,edger,edgeh,edgeb,space_h,space_v);
% pos : [left bottom width height]
[pos]=subplot_pos(nChannels,1,0.01,0.1,0.01,0.1,0,0);
for i = 1:length(pos)
    handles.rnge{ch} = [2*min(cell2mat(data(ch).AP')) 2*max(cell2mat(data(ch).AP'))];
    handles.hAx(i) = subplot(nChannels,1,i,'position',pos{i});
    handles.hSlider(i) = uicontrol('Style','Slider','Units','Normalized','position',[0.9 pos{i}(2) 0.01 pos{i}(4)],'Callback',@sliderCallback,'Value',defThr,'Min',handles.rnge{ch}(1),'Max',handles.rnge{ch}(2)); %'SliderStep',[0.001 0.01],
    handles.thr(i) = defThr;
end
handles.hAx(31) = subplot(nChannels,1,31,'position',pos{31}); % Weird bug

handles.hGlobSlider = uicontrol('Style','Slider','Units','Normalized','Position',[0.95 0.5 0.02 0.05],'Callback',@GlobSliderCallback,'Value',defThr,'Min',handles.rnge{ch}(1),'Max',handles.rnge{ch}(2));
handles.hButtonNext = uicontrol('Style','pushButton','Units','Normalized','Position',[0.97 0.3 0.02 0.05],'Callback',@rangeCallback,'String','>');
handles.hButtonPrev = uicontrol('Style','pushButton','Units','Normalized','Position',[0.93 0.3 0.02 0.05],'Callback',@rangeCallback,'String','<');

handles.sr = sr;
handles.data = data;
handles.lastLoadedSweep = nSweepNeeded;
handles.GlobSlidePrevVal = defThr;
handles.viewRange = [0 10]; % seconds
handles.step = 5; % seconds
handles.totDur = nSweep * sweepDur;

guidata(handles.hF,handles);
replot(handles.hF,1:nChannels);

uiwait;
handles = guidata(handles.hF);
close(handles.hF);


% Full Data loading


for i = 1:nSweep,
    dataS = f32read(fullfile(dataPath,datFiles{i}));
    
    for ch = 1:nChannels,
        dataCh = dataS(ch:nChannels:end); % de-interleave
        LFP = filter(filterLFP,dataCh);
    %         LFP = filtfilt(filterLFP,dataCh);
        AP = filter(filterAP,dataCh);
        
        data(ch).unfiltTrace{i} = dataCh;
        data(ch).LFP{i} = LFP;
        data(ch).AP{i} = AP;
    end
    
    
end




% 
% lpFilt = designfilt('lowpassfir','PassbandFrequency',0.25, ...
%          'StopbandFrequency',0.35,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');
% fvtool(lpFilt)
% dataIn = rand([1000 1]); dataOut = filter(lpFilt,dataIn);



%% trigger spikes

%% get waveforms

%% Save in a sensible format


end

function sliderCallback(hF)

end

function replot(hF,i)

end