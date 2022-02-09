function [center, hausdorff_dists] = buildtboneSSM_propagated()
% function buildtboneSSM(base_dir, part, segment_id, side, index)
% Matlab Code for Atlas Construction
% author: Andy Ding, Alex Lu (Adapted from Runze Han's code)

% create Statistical shape model (SSM) from segmentations with STL 
clear all;
addpath('./CPD');
addpath('./CPD/core');
addpath('./CPD/core/utils');
addpath('./CPD/core/Rigid');
addpath('./CPD/core/Nonrigid');
addpath('./CPD/core/mex');
addpath('./CPD/core/FGT');
addpath('./iso2mesh');
addpath('./HausdorffDist');

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
base_dir = '/Volumes/EXTREME SSD/ANTs-registration';
seg_dir = fullfile(base_dir, 'segmentations');
pred_dir = fullfile(base_dir, 'predictions');
surf_dir = fullfile(base_dir, 'ssm_surfaces');
h5_dir = fullfile(base_dir, 'ssm_H5');

% Create subdirectories if needed
if not(exist(surf_dir, 'dir'))
    fprintf('Creating surfaces directory');
    mkdir(base_dir,'ssm_surfaces');
elseif not(exist(h5_dir, 'dir'))
    fprintf('Creating H5 directory');
    mkdir(base_dir,'ssm_H5');
end

%% 1. Segmentation to Surface Mesh
% Input: segmentations data/Segmentation 1:left tbone, 2:right tbone
% Output: surface mesh vertices: v, faces: f

overwrite_meshes = false;
overwrite_cpd_meshes = false;
mesh_space = 'left-posterior-superior';
include_other_side = true;
viz_cpd_rigid = 0;
viz_cpd_nonrigid = 0;
mutually_exclusive = true;

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
cpd_surf_dir = fullfile(surf_dir, 'cpd');
h5_dir = fullfile(h5_dir, part);
if not(exist(surf_dir, 'dir'))
    fprintf('Creating surfaces directory for part %s\n', part);
    mkdir(surf_dir);
end
if not(exist(h5_dir, 'dir'))
    fprintf('Creating H5 directory for part %s\n', part);
    mkdir(h5_dir);
end
if not(exist(cpd_surf_dir, 'dir'))
    fprintf('Creating CPD surface directory for part %s\n', part);
    mkdir(cpd_surf_dir);
end

template = 153;
[RT, LT] = findSegFileIndices(pred_dir, '.seg.nrrd', side, true, template, include_other_side);
RT(RT == 134 | RT == 145 | RT == 150) = [];
LT(LT == 134 | LT == 150) = [];

if mutually_exclusive
    common_indices = intersect(RT, LT);
    LT = setxor(LT, common_indices);
end

if strcmp(side, 'RT')
    index = RT;
    other_side = 'LT';
    other_side_index = LT;
else
    index = LT;
    other_side = 'RT';
    other_side_index = RT;
end

num_surfaces = length(index);
if include_other_side
    num_surfaces = num_surfaces + length(other_side_index);
end
        
fprintf('Using template %s %d\n', side, template)

for i=0:length(index)
    if i == 0 || index(i) == template
        curr_data = adjustSegPath(seg_dir, '.seg.nrrd', side, template, false, false, false);
        mesh_path = adjustMeshPath(surf_dir, '.vtk', part, side, template, false, false, false, false);
        mesh_zero_path = adjustMeshPath(surf_dir, '.vtk', part, side, template, false, false, true, false);
    else
        curr_data = adjustSegPath(pred_dir, '.seg.nrrd', side, template, true, index(i), false);
        mesh_path = adjustMeshPath(surf_dir, '.vtk', part, side, template, true, index(i), false, false);
        mesh_zero_path = adjustMeshPath(surf_dir, '.vtk', part, side, template, true, index(i), true, false);
    end
    
    if overwrite_meshes || not(exist(mesh_path, 'file')) || not(exist(mesh_zero_path, 'file'))
        if i == 0
            fprintf('Generating %s meshes for %s %d\n', part, side, template)
        else
            fprintf('Generating %s meshes for %s %d\n', part, side, index(i))
        end
        [seg,meta] = nrrdread(curr_data, segment_id); % nrrd path, segment_id
        
        if not(strcmp(meta.space, mesh_space))
            seg = fliplr(flipud(seg));
        end
        
        [Ox,Oy,Oz,sx,sy,sz] = extractMetaFromNrrd(meta); % spacing and orientation of CT
        [nx,ny,nz] = size(seg);
        [v,f] = vol2surf(seg,1:nx,1:ny,1:nz,0.035,1,'simplify'); %mesh simplification ratio, whether to repair mesh holes
        f = f(:,1:3);
        v = [v(:,2) v(:,1) v(:,3)];
        v = [Ox,Oy,Oz]+v.*[sx sy sz]; %convert from matrix indices to physical coordinates
        % write vtks for zero-origined surfaces
        vtkwrite(mesh_zero_path, 'polydata', 'triangle', v(:,1)-Ox, v(:,2)-Oy, v(:,3)-Oz, f, 'Precision', 6);
        % write vtks for surfaces in their original coordinates
        vtkwrite(mesh_path, 'polydata', 'triangle', v(:,1), v(:,2), v(:,3), f, 'Precision', 6);
    else
        if i == 0
            fprintf('%s meshes for %s %d already exist! No need to generate new meshes\n', part, side, template)
        else
            fprintf('%s meshes for %s %d already exist! No need to generate new meshes\n', part, side, index(i))
        end
    end
end

if include_other_side
    for i=1:length(other_side_index)
        curr_data = adjustSegPath(pred_dir, '.seg.nrrd', side, template, true, other_side_index(i), true);
        mesh_path = adjustMeshPath(surf_dir, '.vtk', part, side, template, true, other_side_index(i), false, true);
        mesh_zero_path = adjustMeshPath(surf_dir, '.vtk', part, side, template, true, other_side_index(i), true, true);

        if overwrite_meshes || not(exist(mesh_path, 'file')) || not(exist(mesh_zero_path, 'file'))
            fprintf('Generating flipped %s meshes for %s %d\n', part, other_side, other_side_index(i))
            [seg,meta] = nrrdread(curr_data, segment_id); % nrrd path, segment_id

            if not(strcmp(meta.space, mesh_space))
                seg = fliplr(flipud(seg));
            end
            
            % Flip segments to the same side as template
            seg = fliplr(seg);

            [Ox,Oy,Oz,sx,sy,sz] = extractMetaFromNrrd(meta); % spacing and orientation of CT
            [nx,ny,nz] = size(seg);
            [v,f] = vol2surf(seg,1:nx,1:ny,1:nz,0.035,1,'simplify'); %mesh simplification ratio, whether to repair mesh holes
            f = f(:,1:3);
            v = [v(:,2) v(:,1) v(:,3)];
            v = [Ox,Oy,Oz]+v.*[sx sy sz]; %convert from matrix indices to physical coordinates
            % write vtks for zero-origined surfaces
            vtkwrite(mesh_zero_path, 'polydata', 'triangle', v(:,1)-Ox, v(:,2)-Oy, v(:,3)-Oz, f, 'Precision', 6);
            % write vtks for surfaces in their original coordinates
            vtkwrite(mesh_path, 'polydata', 'triangle', v(:,1), v(:,2), v(:,3), f, 'Precision', 6);
        else
            fprintf('%s meshes for %s %d already exist! No need to generate new meshes\n', part, other_side, other_side_index(i))
            continue
        end
    end
end
% Note: use .vtk for surface mesh. It preserves vertices order. I had
% trouble using .stl.

% Alternatively, one can use external softwares (e.g. MITK)
% load segmentation image -> create smooth polydata -> remeshing

%% 2. Surface Mesh: Obtain Correspondence (CPD)
% Input: segmentation surface mesh [v,f]
% Output: surface mesh with correspondence
%         v_cpd: nv*n each column is the flattened surface vertices

template_mesh = adjustMeshPath(surf_dir, '.vtk', part, side, template, false, false, false, false);
[v0,f0] = vtkread(template_mesh); %the template surface

% Deform the template to each other surface to achieve correspondence
opt_rigid = struct('method','rigid', 'corresp',1, 'normalize',1, 'scale',1, 'viz',viz_cpd_rigid);
opt_nonrigid = struct('method','nonrigid', 'corresp',1, 'normalize',1, 'outliers',0.01, 'viz',viz_cpd_nonrigid, 'beta', 4, 'lambda',3);
v_cpd = cell(num_surfaces,1); %initiate matrix for corresponding surfaces
hausdorff_dists = [];
v_cpd{1} = v0;

% Deform template to surfaces of the same side
for i=2:length(index)
    mesh_cpd_path = adjustMeshPath(cpd_surf_dir, '_cpd.vtk', part, side, template, true, index(i), false, false);
    mesh_path = adjustMeshPath(surf_dir, '.vtk', part, side, template, true, index(i), false, false);
    [vi,~] = vtkread(mesh_path);
    if overwrite_cpd_meshes || not(exist(mesh_cpd_path, 'file'))
        rigid = cpd_register(vi,v0,opt_rigid);
        deform = cpd_register(vi,rigid.Y,opt_nonrigid);
        temp = deform.Y;
        v_cpd{i} = temp;
        vtkwrite(mesh_cpd_path,'polydata','triangle',temp(:,1),temp(:,2),temp(:,3),f0,'Precision',6);
    else
        fprintf('%s CPD meshes for %s %d already exist! No need to generate new meshes\n', part, side, index(i))
        [temp,~] = vtkread(mesh_cpd_path);
        v_cpd{i} = temp;
    end
    hausdorff_dists(i) = HausdorffDist(vi, temp);
end

% Deform template to surfaces of the opposite side that have been flipped
if include_other_side
    for j=1:length(other_side_index)
        mesh_cpd_path = adjustMeshPath(cpd_surf_dir, '_cpd.vtk', part, side, template, true, other_side_index(j), false, true);
        mesh_path = adjustMeshPath(surf_dir, '.vtk', part, side, template, true, other_side_index(j), false, true);
        [vj,~] = vtkread(mesh_path);
        if overwrite_cpd_meshes || not(exist(mesh_cpd_path, 'file'))
            rigid = cpd_register(vj,v0,opt_rigid);
            deform = cpd_register(vj,rigid.Y,opt_nonrigid);
            temp = deform.Y;
            v_cpd{length(index)+j} = temp;
            vtkwrite(mesh_cpd_path,'polydata','triangle',temp(:,1),temp(:,2),temp(:,3),f0,'Precision',6);
        else
            fprintf('%s CPD meshes for %s %d already exist! No need to generate new meshes\n', part, other_side, other_side_index(j))
            [temp,~] = vtkread(mesh_cpd_path);
            v_cpd{length(index)+j} = temp;
        end
        hausdorff_dists(length(index)+j) = HausdorffDist(vj, temp);
    end
end

% Rigidly align each other surface to the template
v_rigid = zeros(numel(v0),num_surfaces);
v_rigid(:,1) = reshape(v0',[numel(v0) 1]); %flatten N*3 surface vertices into a 3N vector
for i=2:num_surfaces
    [~,Z] = procrustes(v0,v_cpd{i});  
    v_rigid(:,i) = reshape(Z',[numel(Z) 1]);
end


%% 3. Build SSM using PCA
[coeff,score,latent] = pca(v_rigid');
pcaMean = mean(v_rigid,2)
mu = reshape(pcaMean,[3 numel(pcaMean)/3])';
center = mean(mu,1);

%% 4. Write SSM to file (.h5)
%pcaMean n*1, coeff n*d, latent d*1, v 3*n, f: 3*nf
h5_path = adjustH5Path(h5_dir, '.h5', part, side, template, false)
writeH5(h5_path,pcaMean,coeff,latent,v0',f0'-1);

%% 5. Visualize final Procrustes alignment
figure;
hold on;
% trimesh(f0, shape(:
scatter3(mu(:,1), mu(:,2), mu(:,3), 'MarkerEdgeColor', 'k', 'DisplayName', 'pcaMean');
for i=1:length(index)
    shape = reshape(v_rigid(:,i), [3 numel(pcaMean)/3])';
%     scatter3(shape(:,1), shape(:,2), shape(:,3), 'MarkerEdgeColor', [rand, rand, rand], 'DisplayName', sprintf('%s %d target', part, index(i)))
    trimesh(f0, shape(:,1), shape(:,2), shape(:,3), 'EdgeColor', 'k', 'FaceColor', [.1 .1 .1], 'FaceAlpha', .05, 'EdgeAlpha', .1, 'DisplayName', sprintf('%s %d target', part, index(i)))
    scatter3(shape(:,1), shape(:,2), shape(:,3), 'MarkerEdgeColor', [rand, rand, rand], 'DisplayName', sprintf('%s %d target', side, index(i)))
end
for i=1:length(other_side_index)
    shape = reshape(v_rigid(:,i), [3 numel(pcaMean)/3])';
%     scatter3(shape(:,1), shape(:,2), shape(:,3), 'MarkerEdgeColor', [rand, rand, rand], 'DisplayName', sprintf('%s %d target', part, index(i)))
    trimesh(f0, shape(:,1), shape(:,2), shape(:,3), 'EdgeColor', 'k', 'FaceColor', [.1 .1 .1], 'FaceAlpha', .05, 'EdgeAlpha', .1, 'DisplayName', sprintf('%s %d target', part, other_side_index(i)))
    scatter3(shape(:,1), shape(:,2), shape(:,3), 'MarkerEdgeColor', [rand, rand, rand], 'DisplayName', sprintf('%s %d target', other_side, other_side_index(i)))
end

legend
hold off;

%% mesh visualization
figure;
trimesh(f0, v0(:,1), v0(:,2), v0(:,3), 'EdgeColor', 'k')

end