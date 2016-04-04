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