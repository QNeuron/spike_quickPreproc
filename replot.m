
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