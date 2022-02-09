function path = adjustH5Path(dir, suffix, part, side, template, zero_origin)
    path = sprintf(fullfile(dir, sprintf('SSM_%s_%s_%d', part, side, template)));
    if zero_origin
        path = strcat(path,'_zero_origin');
    end
    path = strcat(path,suffix);
end