function deform_mesh()
% function deform_mesh()
% Matlab Code for Face Mask Mesh Deformation
% author: Runze Han

% create Statistical shape model (SSM) from segmentations with STL 
addpath('./CPD');
addpath('./CPD/core');
addpath('./CPD/core/utils');
addpath('./CPD/core/Rigid');
addpath('./CPD/core/Nonrigid');
addpath('./CPD/core/mex');
addpath('./CPD/core/FGT');
addpath('./iso2mesh');

%% 0. Housekeeping
base_dir = '/Users/andyding/Documents/Johns Hopkins/OHNS Research/Defacing Head CTs/Datasets';
seg_dir = fullfile(base_dir, 'head_facemasks');
% pred_dir = fullfile(base_dir, 'predictions');
surf_dir = fullfile(base_dir, 'head_facemask_surfaces');

% Create subdirectories if needed
if not(exist(surf_dir, 'dir'))
    fprintf('Creating surfaces directory');
    mkdir(surf_dir)
end

%% 1. Segmentation to Surface Mesh
viz_cpd_rigid = 1;
viz_cpd_nonrigid = 1;
overwrite_meshes = false;
template = 100324;

targets = dir(fullfile(seg_dir, "case-*_THIN_BONE_HEAD_Segmentation.nii.gz"));

for i=1:length(targets)
    target = sscanf(targets(i).name, "case-%d_THIN_BONE_HEAD_Segmentation.nii.gz");
    
    curr_data = fullfile(seg_dir, sprintf("case-%d_THIN_BONE_HEAD_Segmentation.nii.gz", target));
    mesh_path = fullfile(surf_dir, sprintf("case-%d_THIN_BONE_HEAD_Segmentation.vtk", target));
    mesh_zero_path = fullfile(surf_dir, sprintf("case-%d_THIN_BONE_HEAD_Segmentation_zero.vtk", target));
    
    if overwrite_meshes || not(exist(mesh_path, 'file')) || not(exist(mesh_zero_path, 'file'))
        fprintf('Generating mesh for %d\n', target);
        seg = niftiread(curr_data);
        meta = niftiinfo(curr_data);
        [Ox,Oy,Oz,sx,sy,sz] = extractMetaFromNifti(meta); % spacing and orientation of CT
        [nx,ny,nz] = size(seg);
        [v,f] = vol2surf(seg,1:nx,1:ny,1:nz,0.035,1,'simplify'); %mesh simplification ratio, whether to repair mesh holes
        f = f(:,1:3);
        v = [v(:,2) v(:,1) v(:,3)];
        v = transformPointsForward(meta.Transform, v); %convert from matrix indices to physical coordinates
        % write vtks for zero-origined surfaces
        vtkwrite(mesh_zero_path,...
            'polydata', 'triangle', v(:,1)-Ox, v(:,2)-Oy, v(:,3)-Oz, f, 'Precision', 6);
        % write vtks for surfaces in their original coordinates
        vtkwrite(mesh_path,...
            'polydata', 'triangle', v(:,1), v(:,2), v(:,3), f, 'Precision', 6);
    else
        fprintf('Mesh for %d already exists! No need to generate new meshes\n', target);
        % continue
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

num_surfaces = length(targets);
template_mesh = fullfile(surf_dir, sprintf("case-%d_THIN_BONE_HEAD_Segmentation.vtk", template));
[v0,f0] = vtkread(template_mesh); %the template surface

% Deform the template to each other surface to achieve correspondence
opt_rigid = struct('method','rigid', 'corresp',1, 'normalize',1, 'scale',1, 'viz',viz_cpd_rigid);
opt_nonrigid = struct('method','nonrigid', 'corresp',1, 'normalize',1, 'outliers',0.01, 'viz',viz_cpd_nonrigid, 'beta',1, 'lambda',1, 'max_it',40);
v_cpd = cell(num_surfaces,1); %initiate matrix for corresponding surfaces
v_cpd{1} = v0;

% Deform template to surfaces of the same side
for i=1:length(targets)
    target = sscanf(targets(i).name, "case-%d_THIN_BONE_HEAD_Segmentation.nii.gz");
    mesh_cpd_path = fullfile(surf_dir, sprintf("case-%d-%d_THIN_BONE_HEAD_Segmentation.vtk", template, target));
    if template ~= target && ~isfile(mesh_cpd_path)
        fprintf("Registering template %d onto target %d\n", template, target);
        target_mesh = fullfile(surf_dir, sprintf("case-%d_THIN_BONE_HEAD_Segmentation.vtk", target));

        [vi,~] = vtkread(target_mesh);
        rigid = cpd_register(vi,v0,opt_rigid);
        deform = cpd_register(vi,rigid.Y,opt_nonrigid);
        temp = deform.Y;
        v_cpd{i} = temp;
        vtkwrite(mesh_cpd_path,'polydata','triangle',temp(:,1),temp(:,2),temp(:,3),f0,'Precision',6);
    end
end


end