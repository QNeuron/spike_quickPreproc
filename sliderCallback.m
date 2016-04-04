
function sliderCallback(hObj,event)
handles = guidata(hObj);
ch = find(handles.hSlider==hObj);
handles.thr(ch) = get(hObj,'Value');
guidata(hObj,handles);
replotLine(hObj,ch);
end