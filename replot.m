
function replot(hObj,i)
handles = guidata(hObj);

for ch = i
%     axes(handles.hAx(ch));
    y = cell2mat(handles.data(ch).AP');
    t = 0:1/handles.sr:(length(y)-1)/handles.sr;
    rgi = t <= handles.viewRange(2) & t >= handles.viewRange(1);
    t = t(rgi);
    y = y(rgi);
    plot(handles.hAx(ch),t,y);
    hold on;
    handles.hLines(ch) = line([handles.viewRange(1) handles.viewRange(2)],[handles.thr(ch) handles.thr(ch)],'Parent',handles.hAx(ch),'color','red');
    hold off;
    set(handles.hAx(ch),'ytick',[]);
    if ch ~= 32,
        set(handles.hAx(ch),'xtick',[]);
    end
end

guidata(hObj,handles);

end