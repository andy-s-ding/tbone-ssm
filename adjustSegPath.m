function path = adjustSegPath(dir, suffix, side, template, propagation, target, flip)
    if isa(template, 'numeric')
        template = num2str(template);
    end
    if isa(target, 'numeric')
        disp('here')
        target = num2str(target)
    end
    
    if propagation
        path = sprintf(fullfile(dir, sprintf('Segmentation %s %s %s-syn80-demons-downsample300', side, template, target)));
        if flip
            path = strcat(path,'-flipped');
        end
    else
        path = sprintf(fullfile(dir, sprintf('Segmentation %s %s', side, template)));
    end
    path = strcat(path,suffix);
end