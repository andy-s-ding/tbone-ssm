function path = adjustMeshPath(dir, suffix, part, side, template, propagation, target, zero_origin, flip)
    if propagation
        path = sprintf(fullfile(dir, sprintf('%s_%s_%d_%d', part, side, template, target)));
    else
        path = sprintf(fullfile(dir, sprintf('%s_%s_%d', part, side, template)));
    end
    if flip
        path = strcat(path,'_flipped');
    end
    if zero_origin
        path = strcat(path,'_zero_origin');
    end
    path = strcat(path,suffix);
end