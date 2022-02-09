base_dir = '/Volumes/Extreme SSD/Dyna CT Segmentations';
part = 'vestibule_+_cochlea'; % List other segments here
segment_id = 5;
side = 'RT'; %'RT' = right, 'LT' = left

if strcmp(side, 'RT')
    % omit 134 because of dense artifact from cochlear implant
    index = [138, 142, 143, 144, 146, 147, 150, 152, 153]; %test 
else
    index = [138, 143, 144, 145, 146, 147, 148, 150, 151, 152, 169, 171];
end

% right indices = [138, 142, 143, 144, 146, 147, 150, 152, 153] (running list)
% left indices = [138, 143, 144, 145, 146, 147, 148, 150, 151, 152, 169, 171] (running list)

buildtboneSSM(base_dir, part, segment_id, side, index)