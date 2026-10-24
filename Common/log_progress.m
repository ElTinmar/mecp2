function log_progress(varargin)
%%% log_progress(format, ...)
% fprintf to the command window without recording it in the diary
    without_diary(@() fprintf(varargin{:}));
end
