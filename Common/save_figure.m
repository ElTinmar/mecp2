function save_figure(f,basename)
%%% save_figure(f,basename)
% save figure f as basename.pdf and basename.svg, keeping export warnings
% out of the diary
    without_diary(@() print(f,[basename '.pdf'],'-dpdf'));
    without_diary(@() print(f,[basename '.svg'],'-dsvg'));
end
