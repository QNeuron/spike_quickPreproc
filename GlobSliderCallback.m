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