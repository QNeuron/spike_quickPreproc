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