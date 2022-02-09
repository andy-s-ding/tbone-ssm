clear all;
%% Segments
% 1: Bone
% 2: Malleus
% 3: Incus
% 4: Stapes
% 5: Vestibule_and_Cochlea
% 6: Vestibulocochlear_Nerve
% 7: Superior_Vestibular_Nerve
% 8: Inferior_Vestibular_Nerve
% 9: Cochlear_Nerve
% 10: Facial_Nerve
% 11: Chorda_Tympani
% 12: ICA
% 13: Sinus_and_Dura
% 14: Vestibular_Aqueduct
% 15: TMJ
% 16: EAC

seg_names = {
    'Bone',
    'Malleus',
    'Incus',
    'Stapes',
    'Bony Labyrinth',
    'IAC',
    'Superior Vestibular Nerve',
    'Inferior Vestibular Nerve',
    'Cochlear Nerve',
    'Facial Nerve',
    'Chorda',
    'ICA',
    'Sinus + Dura',
    'Vestibular Aqueduct',
    'Mandible',
    'EAC',
};
seg_ids = 1:length(seg_names);
seg_map = containers.Map(seg_ids, seg_names);

%% 0. Housekeeping
groundtruth = false;
name = "Andy";
if strcmp(name, "Andy")
    base_dir = '/Volumes/EXTREME SSD/ANTs-registration';
    image_file_dir = '';
else
    base_dir = '/Users/alex/Documents/Segmentations';
    image_file_dir = '/Users/alex/Documents/images';
end

surf_dir = fullfile(base_dir, 'ssm_surfaces');
h5_dir = fullfile(base_dir, 'ssm_H5');

if groundtruth
    surf_dir = append(surf_dir, "_groundtruth");
    h5_dir = append(h5_dir, "_groundtruth");
end

if not(exist(surf_dir, 'dir'))
    fprintf('Creating surfaces directory');
    mkdir(base_dir,'ssm_surfaces')
elseif not(exist(h5_dir, 'dir'))
    fprintf('Creating H5 directory');
    mkdir(base_dir,'ssm_H5')
end

%% 1. Loading things in
side = 'RT';
segment_id = 5;
if length(segment_id) == 1
    part = seg_map(segment_id);
else
    segment_id = sort(segment_id);
    part = seg_map(segment_id(1));
    for i=2:length(segment_id)
        part = strcat(part, '-', seg_map(segment_id(i)));
    end
end
surf_dir = fullfile(surf_dir, part);
h5_dir = fullfile(h5_dir, part);
pca_surf_dir = fullfile(h5_dir, 'pca_surfaces');

if not(exist(surf_dir, 'dir'))
    fprintf('Creating surfaces directory for part %s', part);
    mkdir(surf_dir);
elseif not(exist(h5_dir, 'dir'))
    fprintf('Creating H5 directory for part %s', part);
    mkdir(h5_dir);
elseif not(exist(pca_surf_dir, 'dir'))
    fprintf('Creating PCA directory for part %s', part);
    mkdir(pca_surf_dir);
end

template = 153;

fSSM = adjustH5Path(h5_dir, '.h5', part, side, template, false)

mu = h5read(fSSM,'/model/mean');
P = h5read(fSSM,'/model/pcaBasis');
faces = h5read(fSSM,'/representer/cells')+1;
mu = reshape(mu,[3 numel(mu)/3])';
P = P';
% size(P)

pcaVariance = h5read(fSSM, '/model/pcaVariance')
pca_var_path = adjustH5Path(h5_dir, '_pca_variance', part, side, template, false);
pca_cum_var_path = adjustH5Path(h5_dir, '_pca_cum_variance', part, side, template, false);

cum_pca_var = ones(length(pcaVariance),1);
for i=1:length(pcaVariance)
    cum_pca_var(i) = sum(pcaVariance(1:i))/sum(pcaVariance);
end

figure;
plot(pcaVariance/sum(pcaVariance), 'marker', 'o');
set(gca, 'XTick', 1:length(pcaVariance));
title(sprintf('Principal Component Analysis: %s', strrep(part, '_', ' and ')))
xlabel('Principal Component Number') 
ylabel('Percentage of Variance Explained')
set(gca,'XTick',0:2:length(pcaVariance)*2);
saveas(gcf, pca_var_path, 'png');
hold off

figure;
plot(cum_pca_var, 'marker', 'o');
set(gca, 'XTick', 1:length(pcaVariance));
title(sprintf('Principal Component Analysis: %s', strrep(part, '_', ' and ')))
xlabel('Principal Component Number') 
ylabel('Cumulative Percentage of Variance Explained')
set(gca,'XTick',0:2:length(pcaVariance)*2);
saveas(gcf, pca_cum_var_path, 'png');
hold off

scalar = 4*cos(linspace(0, 200*pi, 100));
%% janky loop animation 
% f = figure();
% axis vis3d;
% rotate3d on;
% scatter3(mu(:,1), mu(:,2), mu(:,3), 'MarkerEdgeColor', 'k', 'DisplayName', 'pcaMean');
% 
% which_pc = 2;
% 
% for i=1:length(scalar)
%     cla(f)
%     hold on
%     scatter3(mu(:,1), mu(:,2), mu(:,3), 'MarkerEdgeColor', 'k', 'DisplayName', 'pcaMean');
% 
%     PRes = scalar(i)*P(:,which_pc);
%     shape = mu + reshape(PRes,size(mu'))';
%     scatter3(shape(:,1), shape(:,2), shape(:,3), '.', 'MarkerEdgeColor', 'blue', 'DisplayName', sprintf('PC #%d', which_pc), 'MarkerEdgeAlpha', .5);
%     hold off
%     pause(.01)
%     legend   
% end

%% Write out +/- 2 std dev in direction of first 3 principal components + mean

for i=1:3
    for j=-2:2
        shape = mu + reshape(j*sqrt(pcaVariance(i))*P(:,i), size(mu'))';
        mesh_pca_path = adjustH5Path(pca_surf_dir, sprintf('_pca_%d_%d.vtk',i,j), part, side, template, false);
        vtkwrite(mesh_pca_path,'polydata','triangle',shape(:,1), shape(:,2), shape(:,3), faces,'Precision',6)
    end
    
end

%% Write out +/- 2 std dev for first 2 principal components (ensemble) + mean
mode_weights = zeros(size(pcaVariance));
for i=-2:2
    for j=-2:2
        mode_weights(1) = i;
        mode_weights(2) = j;
        shape = mu + reshape(P*(mode_weights.*sqrt(pcaVariance)), size(mu'))';
        mesh_pca_path = adjustH5Path(pca_surf_dir, sprintf('_pca_ensemble_%d_%d.vtk',i,j), part, side, template, false);
        vtkwrite(mesh_pca_path,'polydata','triangle',shape(:,1), shape(:,2), shape(:,3), faces,'Precision',6)
    end
    
end


%% Write out the mean structure

mesh_pca_mean_path = adjustH5Path(pca_surf_dir, '_pca_mean.vtk', part, side, template, false);
vtkwrite(mesh_pca_mean_path,'polydata','triangle',mu(:,1), mu(:,2), mu(:,3), faces,'Precision',6)


%% interactive slider visualization
% adapted from https://www.mathworks.com/help/control/ug/build-app-with-interactive-plot-updates.html

f = figure;
hold on;
scatter3(mu(:,1), mu(:,2), mu(:,3), 'MarkerEdgeColor', 'k', 'DisplayName', 'pcaMean');

pc_coeff = 0;
i = 1;
shape = mu + reshape(pc_coeff*sqrt(pcaVariance(i))*P(:,i), size(mu'))';
h=scatter3(shape(:,1), shape(:,2), shape(:,3), 'MarkerEdgeColor', 'b', 'DisplayName', 'in PC space');
hold off

legend

b1 = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',pc_coeff, 'min',-3, 'max',3,...
              'String', 'coeff pc1',...
              'Tag','slider1',...
              'UserData',struct('coeff',0),...
              'Callback', {@scatter_callback, mu, P, pcaVariance});
          
b2 = uicontrol('Parent',f,'Style','slider','Position',[81,20,419,23],...
              'value',pc_coeff, 'min',-3, 'max',3,...
              'String', 'coeff pc2',...
              'Tag','slider2',...
              'UserData',struct('coeff',0),...
              'Callback', {@scatter_callback, mu, P, pcaVariance});
          
b3 = uicontrol('Parent',f,'Style','slider','Position',[81,0,419,23],...
              'value',pc_coeff, 'min',-3, 'max',3,...
              'String', 'coeff pc3',...
              'Tag','slider3',...
              'UserData',struct('coeff',0),...
              'Callback', {@scatter_callback, mu, P, pcaVariance});
          
          
% bgcolor = f.Color;
% bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
%                 'String','-3','BackgroundColor',bgcolor);
% bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
%                 'String','3','BackgroundColor',bgcolor);
% bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
%                 'String','Coeff','BackgroundColor',bgcolor);

% b.Callback = @(es,ed) updateSystem(h,tf(wn^2,[1,2*(es.Value)*wn,wn^2])); 

function scatter_callback(src, event, mu, P, pcaVariance)
    switch src.String 
        case 'coeff pc1'
            i = 1;
        case 'coeff pc2'
            i = 2;
        case 'coeff pc3'
            i = 3;
    end
    cla
    hold on;
    scatter3(mu(:,1), mu(:,2), mu(:,3), 'MarkerEdgeColor', 'k', 'DisplayName', 'pcaMean');
    
    % assemble the coefficient vector by determining the global state of
    % all sliders
    b = zeros(3, 1);
    for j=1:3
        if j == i
            continue
        end
        h = findobj('Tag',sprintf('slider%d',j));
        data = h.UserData;
        b(j) = h.Value;
    end
    b(i) = src.Value;
    size(b);
    size(pcaVariance);
    shape = mu + reshape(P(:,1:3)*(b.*sqrt(pcaVariance(1:3))), size(mu'))';
    scatter3(shape(:,1), shape(:,2), shape(:,3), 'MarkerEdgeColor', 'b', 'DisplayName', 'in PC space');
    title(sprintf('$$%.2f * \\sqrt(\\lambda_1) * PC1 + %.2f * \\sqrt(\\lambda_2) * PC2 + %.2f * \\sqrt(\\lambda_3) * PC3 + \\mu$$', b(1), b(2), b(3)), 'interpreter', 'latex');
    hold off;
    src.UserData = struct('coeff', src.Value);

end