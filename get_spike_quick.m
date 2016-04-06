function [spikeMat, stimNames, waveforms, APTraces, LFPTraces, thresholds ] = get_spike_quick(dataPath,thresholds)
% spikeMat = [absolute time, relative time, stimulus #, channel]


%% Get data ID

if ~exist('dataPath','var') || isempty(dataPath),
    dataPath = cd;
end
if ~exist('thresholds','var'),
    thresholds = [];
end

files = excludeDots(dir(dataPath));
datFiles = {files(arrayfun(@(x)(~isempty(strfind(x.name,'.f32.dat'))),files)).name};


load(fullfile(dataPath,'gridInfo.mat'))
nChannels = expt.nChannels;
nSweep = grid.nSweepsDesired;
sr = expt.dataDeviceSampleRate;
stimNames = grid.stimFiles;

%% low pass / band pass filter definition
lfpFreqLim = 500; % Hz

filterLFP = designfilt('lowpassiir', ...        % Response type
       'PassbandFrequency',lfpFreqLim, ...     % Frequency constraints
       'StopbandFrequency',lfpFreqLim+150, ...
       'PassbandRipple',4, ...          % Magnitude constraints
       'StopbandAttenuation',65, ...
       'DesignMethod','butter', ...      % Design method
       'MatchExactly','passband', ...   % Design method options
       'SampleRate',sr);               % Sample rate

filterAP = designfilt('bandpassiir', ...       % Response type
       'StopbandFrequency1',lfpFreqLim, ...    % Frequency constraints
       'PassbandFrequency1',lfpFreqLim+150, ...
       'PassbandFrequency2',10000, ...
       'StopbandFrequency2',11000, ...
       'StopbandAttenuation1',65, ...   % Magnitude constraints
       'PassbandRipple',1, ...
       'StopbandAttenuation2',65, ...
       'DesignMethod','butter', ...      % Design method
       'MatchExactly','passband', ...   % Design method options
       'SampleRate',sr);               % Sample rate
% fvtool(filterLFP);


%% Data loading
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
sweepLength = length(dataS)/nChannels; % points
% nSweepNeeded = ceil(plotDur/sweepDur);

for i = 2:nSweep %nSweepNeeded, % Load more seep if more are needed
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

%% GUI block

if isempty(thresholds),
    
    defThr = 0;
    handles.hF = figure('units','normalized','outerposition',[0 0 1 1]);
    % [pos]=subplot_pos(nChannels,1,edgel,edger,edgeh,edgeb,space_h,space_v);
    % pos : [left bottom width height]
    [pos]=subplot_pos(nChannels,1,0.01,0.1,0.01,0.1,0,0);
    for i = 1:length(pos)
        handles.rnge{i} = [2*min(cell2mat(data(i).AP')) 2*max(cell2mat(data(i).AP'))];
        %     handles.ylim{ch} = [1.5*min(cell2mat(data(ch).AP')) 1.5*max(cell2mat(data(ch).AP'))];
        handles.hAx(i) = subplot(nChannels,1,i,'position',pos{i});
        handles.hSlider(i) = uicontrol('Style','Slider','Units','Normalized','position',[0.9 pos{i}(2) 0.01 pos{i}(4)],'Callback',@sliderCallback,'Value',defThr,'Min',handles.rnge{i}(1),'Max',handles.rnge{i}(2)); %'SliderStep',[0.001 0.01],
        handles.thr(i) = defThr;
    end
    handles.hAx(31) = subplot(nChannels,1,31,'position',pos{31}); % Weird bug
    handles.ylim = [min(cell2mat(handles.rnge)) max(cell2mat(handles.rnge))];
    
    handles.hGlobSlider = uicontrol('Style','Slider','Units','Normalized','Position',[0.95 0.5 0.02 0.05],'Callback',@GlobSliderCallback,'Value',defThr,'Min',handles.rnge{1}(1),'Max',handles.rnge{1}(2));
    handles.hButtonNext = uicontrol('Style','pushButton','Units','Normalized','Position',[0.97 0.3 0.02 0.05],'Callback',@rangeCallback,'String','>');
    handles.hButtonPrev = uicontrol('Style','pushButton','Units','Normalized','Position',[0.93 0.3 0.02 0.05],'Callback',@rangeCallback,'String','<');
    handles.hButtonZIn = uicontrol('Style','pushButton','Units','Normalized','Position',[0.93 0.4 0.02 0.05],'Callback',@zoomCallback,'String','Z in');
    handles.hButtonZout = uicontrol('Style','pushButton','Units','Normalized','Position',[0.97 0.4 0.02 0.05],'Callback',@zoomCallback,'String','Z out');
    handles.hButtonDone = uicontrol('Style','pushButton','Units','Normalized','Position',[0.95 0.2 0.02 0.05],'Callback',@(x,y)(uiresume),'String','Done');
    
    handles.sr = sr;
    handles.data = data;
    % handles.lastLoadedSweep = nSweepNeeded;
    % handles.nSweepNeeded = nSweepNeeded;
    handles.GlobSlidePrevVal = defThr;
    handles.viewRange = [0 5]; % seconds
    handles.step = 3; % seconds
    handles.totDur = nSweep * sweepDur;
    handles.nSweep = nSweep;
    handles.sweepLength = sweepLength;
    
    guidata(handles.hF,handles);
    replot(handles.hF,1:nChannels);
    
    uiwait;
    handles = guidata(handles.hF); % Gets the last version of handles
    close(handles.hF);
    
    thresholds = handles.thr;
    
end
%% Spike extraction
mVthreshold = 0.005; % V (?)
durThreshold = 15 * 10^-3; % msec
waveFormDur = 10 * 10^-3; % msec
waveFormBin = ceil(waveFormDur*sr);

c = 1;
for ch = 1:nChannels,
    for s = 1:nSweep,
        % trigger spikes
        [ bin_peak_max val_peak_max peak_val peak_abs ] = peak_detector_general(data(ch).AP{s}, 'thr', thresholds(ch));
        
        % Discard obvious shit (to long & to strong dV)
        GoodInd = (val_peak_max < mVthreshold) & (diff(peak_abs')/sr < durThreshold);
        bin_peak_max = bin_peak_max(GoodInd);
        val_peak_max = val_peak_max(GoodInd);
        peak_val = peak_val(GoodInd);
        peak_abs = peak_abs(GoodInd,:);
        
%         spikeMat, stimNames, waveforms,
        for i = 1:length(bin_peak_max),
            spikeMat(c,1:4) = [(bin_peak_max(i)+(s-1)*sweepLength)/sr bin_peak_max(i)/sr grid.randomisedGrid(s) ch];
            ind = bin_peak_max(i)-floor(waveFormBin/2):bin_peak_max(i)+ceil(waveFormBin/2)-1;
            waveforms(c,1:waveFormBin) = data(ch).AP{s}(ind);
            c = c+1;
        end
        % get waveforms
        
        
        
    
    end
    APTraces(ch,1:nSweep) = data(ch).AP;
    LFPTraces(ch,1:nSweep) = data(ch).LFP;
end





% Full Data loading
% 
% for i = 1:nSweep,
%     dataS = f32read(fullfile(dataPath,datFiles{i}));
%     
%     for ch = 1:nChannels,
%         dataCh = dataS(ch:nChannels:end); % de-interleave
%         LFP = filter(filterLFP,dataCh);
%     %         LFP = filtfilt(filterLFP,dataCh);
%         AP = filter(filterAP,dataCh);
%         
%         data(ch).unfiltTrace{i} = dataCh;
%         data(ch).LFP{i} = LFP;
%         data(ch).AP{i} = AP;
%     end
%     
%     
% end




% % 
% lpFilt = designfilt('lowpassfir','PassbandFrequency',0.25, ...
%          'StopbandFrequency',0.35,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');
% fvtool(lpFilt)
% dataIn = rand([1000 1]); dataOut = filter(lpFilt,dataIn);



end

function replot(hObj,i)
handles = guidata(hObj);

steps = 1:handles.sweepLength:handles.sweepLength*(handles.nSweep+1);
wantedPos = (handles.viewRange * handles.sr)+1;
wantedLim = [find(steps <= wantedPos(1),1,'first') find(steps > wantedPos(2)-1,1,'first')-1];

for ch = i
%     axes(handles.hAx(ch));
    y = cell2mat(handles.data(ch).AP(wantedLim(1):wantedLim(2))');
%     t = 0:1/handles.sr:(length(y)-1)/handles.sr;
    t = ((wantedLim(1)-1)*handles.sweepLength+1)/handles.sr:1/handles.sr:(wantedLim(2)*handles.sweepLength)/handles.sr;
    rgi = t <= handles.viewRange(2) & t >= handles.viewRange(1);
    t = t(rgi);
    y = y(rgi);
    plot(handles.hAx(ch),t,y);
    hold on;
    handles.hLines(ch) = line([handles.viewRange(1) handles.viewRange(2)],[handles.thr(ch) handles.thr(ch)],'Parent',handles.hAx(ch),'color','red');
    hold off;
    set(handles.hAx(ch),'ytick',[]);
    set(handles.hAx(ch),'YLim',handles.ylim);
    if ch ~= 32,
        set(handles.hAx(ch),'xtick',[]);
    end
end

guidata(hObj,handles);

end

function zoomCallback(hObj,event)

handles = guidata(hObj);

switch get(hObj,'String')
    case 'Z in'
%         if handles.viewRange(2) + handles.step > handles.totDur,
%             return;
%         end
        handles.ylim = handles.ylim - 0.2 * handles.ylim;
    case 'Z out'

        handles.ylim = handles.ylim + 0.2 * handles.ylim;
end

for ch = 1:length(handles.hAx),
    set(handles.hAx(ch),'YLim',handles.ylim)
    if handles.thr(ch) > handles.ylim(1) && handles.thr(ch) < handles.ylim(2)
        set(handles.hSlider(ch),'Min',handles.ylim(1),'Max',handles.ylim(2));
        set(handles.hGlobSlider,'Min',handles.ylim(1),'Max',handles.ylim(2));
    end
end

guidata(handles.hF,handles);

end

function rangeCallback(hObject,event)
handles = guidata(hObject);

switch get(hObject,'String')
    case '>'
        if handles.viewRange(2) + handles.step > handles.totDur,
            return;
        end
        handles.viewRange = handles.viewRange + handles.step;
    case '<'
        if handles.viewRange(1) - handles.step < 0,
            return;
        end
        handles.viewRange = handles.viewRange - handles.step;
end
guidata(handles.hF,handles);
replot(handles.hF,1:length(handles.hAx));

end

function GlobSliderCallback(hObj,event)
handles = guidata(hObj);
step = get(hObj,'Value') - handles.GlobSlidePrevVal;
handles.GlobSlidePrevVal = get(hObj,'Value');
for i = 1:length(handles.hAx)
    % slider value
    set(handles.hSlider(i),'Value',handles.thr(i) + step); 
end
% handles.thr
handles.thr = handles.thr + step;

guidata(handles.hF,handles);

% replot lines
replotLine(handles.hF,1:length(handles.hAx));


end

function replotLine(hObj,i)
handles = guidata(hObj);

for ch = i
    set(handles.hLines(ch),'YData',[handles.thr(ch) handles.thr(ch)]);
end

% guidata(hObj,handles);

end

function sliderCallback(hObj,event)
handles = guidata(hObj);
ch = find(handles.hSlider==hObj);
handles.thr(ch) = get(hObj,'Value');
guidata(hObj,handles);
replotLine(hObj,ch);
end