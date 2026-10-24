function cleanup = start_table(filename)
%%% cleanup = start_table(filename)
% record command window output to filename (overwritten) until cleanup is
% cleared, i.e. when the calling function returns or errors
    if isfile(filename), delete(filename); end
    diary(filename);
    cleanup = onCleanup(@() diary('off'));
end
