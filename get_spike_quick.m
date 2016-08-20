function [spikeMat, stimNames, waveforms, APTraces, LFPTraces, thresholds, expInfo ] = get_spike_quick(varargin)
% [spikeMat, stimNames, waveforms, APTraces, LFPTraces, thresholds ] =
% get_spike_quick(dataPath,thresholds,Name-value pair arguments)
% get_spike_quick allows to filter the raw data from Benware, to
% manually set a spike threshold and to extract the corresponding spikes.
%
% inputs :  All inputs are optionnal.
%           - dataPath : directory containing the raw data. DEF : cd
%           - thresholds : threshold values, if already known. This will
%           prevent the GUI to pop up and apply the threshold to the data.
%
% Name-value pair arguments :
%           - lfpFreqLim : cutoff frequency to separate LFP from MUA (in Hz)
%           - plotDur : Duration of the data you want to analyze (in
%           seconds). Empty value will apply the function to the entire
%           dataset. Default: 5 sec
%           - waveFormDur : trigged spike waveform duration (in msec). Default 10
%           msec, centered on the local maximal trigged value.
%           - save : save extracted spikes. Def: false
%
% Outputs :
%           - spikeMat : Trigged spike information. [absolute time, relative
%           time, stimulus #, sweep #, channel]
%           - stimNames : Sweep names, taken from the grid info
%           - waveforms : waveforms of the trigged spikes.
%           - APTraces : complete MUA trace (channel x sweep)
%           - LFPTraces : complete LFP trace (channel x sweep)
%           - thresholds : spike thresholds selected by the user
%           - expInfo : struct containing grid and expt from Benware
%

%% Parse inputs
p = inputParser;
% Default values
def_dataPath = cd;
def_thresholds = [];
def_lfpFreqLim = 500; % Hz
def_plotDur = 5; % seconds
def_waveFormDur = 10 * 10^-3; % msec
def_save = false;

addOptional(p,'dataPath',def_dataPath,@isstr);
addOptional(p,'thresholds',def_thresholds,@isnumeric);
addParameter(p,'lfpFreqLim',def_lfpFreqLim,@isnumeric);
addParameter(p,'plotDur',def_plotDur,@isnumeric);
addParameter(p,'waveFormDur',def_waveFormDur,@isnumeric);
addParameter(p,'save',def_save,@islogical);

parse(p,varargin{:}); % Parse inputs

opt = struct;
for i = 1:length(p.Parameters),
    opt.(p.Parameters{i}) = p.Results.(p.Parameters{i});
end


%% Get data ID
%
% if ~exist('dataPath','var') || isempty(dataPath),
%     dataPath = cd;
% end
% if ~exist('thresholds','var'),
%     thresholds = [];
% end

files = excludeDots(dir(opt.dataPath));
datFiles = {files(arrayfun(@(x)(~isempty(strfind(x.name,'.f32.dat'))),files)).name};

grid = struct; expt = struct;
load(fullfile(opt.dataPath,'gridInfo.mat'))
nChannels = expt.nChannels;
nSweep = grid.nSweepsDesired;
sr = expt.dataDeviceSampleRate;
stimNames = grid.stimFiles;
expInfo.grid = grid;
expInfo.exp = expt;
expInfo.dataPath = opt.dataPath;

%% low pass / band pass filter definition

filterLFP = designfilt('lowpassiir', ...        % Response type
    'PassbandFrequency',opt.lfpFreqLim, ...     % Frequency constraints
    'StopbandFrequency',opt.lfpFreqLim+150, ...
    'PassbandRipple',4, ...          % Magnitude constraints
    'StopbandAttenuation',65, ...
    'DesignMethod','butter', ...      % Design method
    'MatchExactly','passband', ...   % Design method options
    'SampleRate',sr);               % Sample rate

filterAP = designfilt('bandpassiir', ...       % Response type
    'StopbandFrequency1',opt.lfpFreqLim, ...    % Frequency constraints
    'PassbandFrequency1',opt.lfpFreqLim+150, ...
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

fprintf('Data loading...\n');

dataS = f32read(fullfile(opt.dataPath,datFiles{1})); % Load first sweep to get sweep length
for ch = 1:nChannels,
    dataCh = dataS(ch:nChannels:end); % de-interleave
    LFP = filter(filterLFP,dataCh);
    AP = filter(filterAP,dataCh);
    
    data(ch).unfiltTrace{1} = dataCh;
    data(ch).LFP{1} = LFP;
    data(ch).AP{1} = AP;
end
sweepDur = length(dataS)/nChannels/sr; % seconds
sweepLength = length(dataS)/nChannels; % points
expInfo.sweepLength = sweepLength;

if ~isempty(opt.plotDur)
    nSweepNeeded = ceil(opt.plotDur/sweepDur);
    for i = 2:nSweepNeeded % Load a subset of the data
        dataS = f32read(fullfile(opt.dataPath,datFiles{i}));
        
        for ch = 1:nChannels,
            dataCh = dataS(ch:nChannels:end); % de-interleave
            LFP = filter(filterLFP,dataCh);
            AP = filter(filterAP,dataCh);
            
            data(ch).unfiltTrace{i} = dataCh;
            data(ch).LFP{i} = LFP;
            data(ch).AP{i} = AP;
        end
        
        
    end
    
else
    for i = 2:nSweep % load all the data
        dataS = f32read(fullfile(opt.dataPath,datFiles{i}));
        
        for ch = 1:nChannels,
            dataCh = dataS(ch:nChannels:end); % de-interleave
            LFP = filter(filterLFP,dataCh);
            AP = filter(filterAP,dataCh);
            
            data(ch).unfiltTrace{i} = dataCh;
            data(ch).LFP{i} = LFP;
            data(ch).AP{i} = AP;
        end
    end 
end

if ~isempty(i),
    nSweepLoaded = i;
else nSweepLoaded = 1;
end
expInfo.nSweepLoaded = nSweepLoaded;

%% GUI block

if isempty(opt.thresholds),
    
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
        handles.hChkBox(i) = uicontrol('Style','CheckBox','Units','Normalized','Position',[0.91 pos{i}(2) 0.01 pos{i}(4)],'Callback',@chkBoxCallback,'Value',1);
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
    handles.SelectedChannels = 1:nChannels;
    
    replot(handles.hF,1:nChannels);
    
    uiwait; % Wait for the "done" button to be pressed
    close(handles.hF);
    
    thresholds = handles.thr;
    thresholds(thresholds==0) = NaN;
    SelCh = handles.SelectedChannels;
else
    thresholds = opt.thresholds;
    SelCh = find(~isnan(thresholds));
end

%% Spike extraction
fprintf('Extracting Spikes...\n');
mVthreshold = 0.005; % V (?)
durThreshold = 15 * 10^-3; % msec
waveFormBin = ceil(opt.waveFormDur*sr);
abscissa = (0:waveFormBin-1)/sr;

c = 1;
spikeMat = [];
waveforms = [];
TspikeMat = zeros(5000,4);
Twaveforms = zeros(5000,waveFormBin);
APTraces = zeros(nChannels,nSweepLoaded*sweepLength);
LFPTraces = zeros(nChannels,nSweepLoaded*sweepLength);
for ch = SelCh,
% for ch = 1:nChannels,
    for s = 1:length(data(1).AP), % Extract AP from the loaded set
        % trigger spikes
        [ bin_peak_max, val_peak_max, ~, peak_abs ] = peak_detector_general(data(ch).AP{s}, 'thr', thresholds(ch));

        % Discard obvious shit (too long | too strong dV | too early | too late)
        GoodInd = (val_peak_max < mVthreshold) & (diff(peak_abs')/sr < durThreshold)' & (bin_peak_max > waveFormBin/2+1) & (bin_peak_max < (sweepLength-waveFormBin/2-1));
        bin_peak_max = bin_peak_max(GoodInd);
%         val_peak_max = val_peak_max(GoodInd);
%         peak_val = peak_val(GoodInd);
%         peak_abs = peak_abs(GoodInd,:);
        
        %         spikeMat, stimNames, waveforms,
        for i = 1:length(bin_peak_max),
            TspikeMat(c,1:5) = [(bin_peak_max(i)+(s-1)*sweepLength)/sr bin_peak_max(i)/sr grid.randomisedGrid(s) s ch];
            ind = bin_peak_max(i)-floor(waveFormBin/2):bin_peak_max(i)+ceil(waveFormBin/2)-1;
            Twaveforms(c,1:waveFormBin) = data(ch).AP{s}(ind);
            c = c+1;
            
            if c == 5001,
                c = 1;
                spikeMat = [spikeMat; TspikeMat];
                waveforms = [waveforms; Twaveforms];
                
            end
            
        end
        % get waveforms
        
    end
    APTraces(ch,1:nSweepLoaded*sweepLength) = cell2mat(cellfun(@(x)(x'),data(ch).AP,'UniformOutput',false));
    LFPTraces(ch,1:nSweepLoaded*sweepLength) = cell2mat(cellfun(@(x)(x'),data(ch).LFP,'UniformOutput',false));
end

% cut the zeros out
spikeMat = spikeMat(1:end-(5000-c),:);
waveforms = waveforms(1:end-(5000-c),:);

% Saving spikes
if opt.save,
    fprintf('Saving...\n');
    SpikeMatNames = {'absolute time', 'relative time', 'stimulus #', 'sweep #', 'channel'};
    save(fullfile(opt.dataPath,sprintf('spikes_%s_P%d_N%d_%s',expInfo.exp.userName,expInfo.exp.penetrationNum,expInfo.exp.exptNum,expInfo.grid.name)) ...
        ,'spikeMat', 'stimNames', 'thresholds', 'expInfo','SpikeMatNames');
    save(fullfile(opt.dataPath,sprintf('spikes_waveforms_%s_P%d_N%d_%s',expInfo.exp.userName,expInfo.exp.penetrationNum,expInfo.exp.exptNum,expInfo.grid.name)) ...
        ,'waveforms', 'expInfo', 'abscissa');
end
fprintf('Done.\n');

%% GUI Subfunctions

    function replot(hObj,i)
        
        steps = 1:handles.sweepLength:handles.sweepLength*(handles.nSweep+1);
        wantedPos = (handles.viewRange * handles.sr)+1;
        wantedLim = [find(steps <= wantedPos(1),1,'first') find(steps > wantedPos(2)-1,1,'first')-1];
        
%         for ch = handles.SelectedChannels
        for ch = i
            y = cell2mat(handles.data(ch).AP(wantedLim(1):wantedLim(2))');
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
        
    end


    function zoomCallback(hObj,event)
        
        switch get(hObj,'String')
            case 'Z in'
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
        
    end

    function rangeCallback(hObject,event)
        
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
        replot(handles.hF,1:length(handles.hAx));
        
    end

    function GlobSliderCallback(hObj,event)
        
        step = get(hObj,'Value') - handles.GlobSlidePrevVal;
        handles.GlobSlidePrevVal = get(hObj,'Value');
        for i = 1:length(handles.hAx)
            % slider value
            set(handles.hSlider(i),'Value',handles.thr(i) + step);
        end
        handles.thr = handles.thr + step;
        
        % replot lines
        replotLine(handles.hF,1:length(handles.hAx));
        
        
    end

    function replotLine(hObj,i)
        
        for ch = i
            set(handles.hLines(ch),'YData',[handles.thr(ch) handles.thr(ch)]);
        end
        
    end

    function sliderCallback(hObj,event)
        
        ch = find(handles.hSlider==hObj);
        handles.thr(ch) = get(hObj,'Value');
        replotLine(hObj,ch);
    end
    
    function chkBoxCallback(hObj,event)
        ch = find(handles.hChkBox==hObj);
        switch get(hObj,'Value'),
            case 1
                handles.SelectedChannels = sort([handles.SelectedChannels ch]);
            case 0
                handles.SelectedChannels = handles.SelectedChannels(handles.SelectedChannels~=ch);
        end
    end
end
