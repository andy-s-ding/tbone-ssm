function [side_indices, other_side_indices] = findSegFileIndices(seg_dir, suffix, side, propagation, template, include_other_side) 
    if side == 'RT'
        other_side = 'LT';
    else
        other_side = 'RT';
    end
    if propagation
        adjustSegPath(seg_dir, suffix, side, template, propagation, '*', false)
        side_files = dir(adjustSegPath(seg_dir, suffix, side, template, propagation, '*', false));
        side_names = {side_files.name};
        side_indices = cellfun(@(x) extractPropagatedIndex(x, propagation), side_names, 'UniformOutput', false);
        if include_other_side
            other_side_files = dir(adjustSegPath(seg_dir, suffix, side, template, propagation, '*', true));
            other_side_names = {other_side_files.name};
            other_side_indices = cellfun(@(x) extractPropagatedIndex(x, propagation), other_side_names, 'UniformOutput', false);
        else
            other_side_indices = {};
        end
    else
        adjustSegPath(seg_dir, suffix, side, template, propagation, false, false)    
        side_files = dir(adjustSegPath(seg_dir, suffix, side, '*', propagation, false, false));
        side_names = {side_files.name};
        side_indices = cellfun(@(x) extractPropagatedIndex(x, propagation), side_names, 'UniformOutput', false);
        if include_other_side
            other_side_files = dir(adjustSegPath(seg_dir, suffix, other_side, '*', propagation, false, false));
            other_side_names = {other_side_files.name};
            other_side_indices = cellfun(@(x) extractPropagatedIndex(x, propagation), other_side_names, 'UniformOutput', false);
        else
            other_side_indices = {};
        end
    end
    side_indices = cell2mat(side_indices);
    other_side_indices = cell2mat(other_side_indices);
end

function index = extractPropagatedIndex(filename, propagation)
   split_name = split(filename, ' ');
   target_phrase = split_name{end};
   if propagation
       target_split = split(target_phrase, '-');
   else
       target_split = split(target_phrase, '.');
   end
   index = str2num(target_split{1});
end