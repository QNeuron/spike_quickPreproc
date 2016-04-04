function replotLine(hObj,i)
handles = guidata(hObj);

for ch = i
    set(handles.hLines(ch),'YData',[handles.thr(ch) handles.thr(ch)]);
end

% guidata(hObj,handles);

end