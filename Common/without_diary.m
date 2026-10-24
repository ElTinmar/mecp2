function without_diary(fun)
%%% without_diary(fun)
% run fun without recording its command window output in the diary
    recording = strcmp(get(0,'Diary'),'on');
    diary off
    fun();
    if recording, diary on; end
end
