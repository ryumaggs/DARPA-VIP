function matlab_to_hdf5_write(GMRA, imgOpts, filename)

% Use the high-level interface to write Matlab GMRA data out to an HDF5 
% file for import into web-based IPCA visualization. Low-level interface
% used for file reset and hard link construction

% Going to filter out any nodes which have less than this number of points
NODE_PTS_LOWER_LIMIT = 10;

fprintf('(re)Creating HDF5 file, %s\n', filename);

% Recreate file each time
fcpl = H5P.create('H5P_FILE_CREATE');
fapl = H5P.create('H5P_FILE_ACCESS');
fid = H5F.create(filename, 'H5F_ACC_TRUNC', fcpl, fapl);
H5F.close(fid);

%% Original data

% HACK: The MNIST data needs to be rotated and flipped to show up correctly
% in the GUIs, so I'm doing this fragile test that if the original data is
% 784 long (28x28) it must be MNIST...
if (size(GMRA.X, 1) == 784)
    GMRA.X = MNIST_image_flip(GMRA.X);
end

fprintf('Writing original data\n');

% Original data : want [D x N]
h5create(filename, '/original_data/data', size(GMRA.X'));
h5write(filename, '/original_data/data', GMRA.X');

% Data type and attributes
% TODO: Need to figure out a good way to specify options
if imgOpts.imageData,
    h5writeatt(filename, '/original_data', 'dataset_type', 'image');
    h5writeatt(filename, '/original_data', 'image_n_rows', int64(imgOpts.imR));
    h5writeatt(filename, '/original_data', 'image_n_columns', int64(imgOpts.imC));
end

%% Per-point original data Labels

% imgOpts will probalby have to evolve into a more general metadata store
%   about the original data

% imgOpts.Labels should be an [n_points, n_labels] array of label data
%        .LabelNames ? a cell array of strings {n_labels}
%           NOTE: for now names must be unique and not contain spaces because
%           they'll be used as the name of the hdf5 path!!
%           TODO: Move to unique ID as path and explicit Name field...
%        .LabelDataTypes ? a cell array of strings with the data types
%           (e.g. 'int32', 'double'...)
%        .LabelVariableTypes ? a cell array of strings listing
%           whether the variable is 'categorical' or 'continuous'
%        .LabelDescriptions (optional) ? descriptions of what the label is
%           or how it was generated
%        .LabelKeys (not used yet) ? a mapping between categorical integer label
%           values and more human-readable strings (not totally worked out
%           yet how to store...)

fprintf('Writing original data labels\n');

for ii = 1:size(imgOpts.Labels, 2)
    label_data = imgOpts.Labels(:,ii);
    label_name = imgOpts.LabelNames{ii};
    label_data_type = imgOpts.LabelDataTypes{ii};
    label_variable_type = imgOpts.LabelVariableTypes{ii};
    label_description = label_name;
    if isfield(imgOpts, 'LabelDescriptions'),
        label_description = imgOpts.LabelDescriptions{ii};
    end
    
    label_path = ['/original_data/labels/' label_name];
    
    h5create(filename, label_path, length(label_data), 'Datatype', label_data_type);
    % not totally sure I shouldn't be using typecast()...
    h5write(filename, label_path, cast(label_data, label_data_type));
    h5writeatt(filename, label_path, 'description', label_description);
    h5writeatt(filename, label_path, 'variable_type', label_variable_type);
    % h5writeatt(filename, label_path, 'key', {{'1', 'digit1'}, {'2', 'digit2'}});
end

%% Per-point original data Diffusion graph (only eigenvectors for now...)

if isfield(GMRA, 'Graph') && isfield(GMRA.Graph, 'EigenVecs'),
    fprintf('Writing original data diffusion embedding\n');
    
    eigenvecs = GMRA.Graph.EigenVecs;
    eigenvals = GMRA.Graph.EigenVals;
    h5create(filename, '/original_data/diffusion_graph/eigenvectors', size(eigenvecs), 'Datatype', 'double');
    h5write(filename, '/original_data/diffusion_graph/eigenvectors', eigenvecs);
    h5create(filename, '/original_data/diffusion_graph/eigenvalues', size(eigenvals), 'Datatype', 'double');
    h5write(filename, '/original_data/diffusion_graph/eigenvalues', eigenvals);
    
    % TODO: Copy over all GMRA.Graph.Opts options to attributes of /original_data/diffusion_graph
end

%% TODO: Need Per-node Labels and Diffusion Embedding section


%% Full tree by ID

% Excluding nodes with n_pts lower than NODE_PTS_LOWER_LIMIT.
% Will be easier to do children if keep list of excluded nodes
% NOTE: zero-based!
excluded_nodes = [];

fprintf('Writing nodes by ID\n');

% First pass, store nodes by ID
n_total_nodes = length(GMRA.cp);
nodes_root = '/full_tree/nodes/';
for ii = 1:n_total_nodes,
    
    if GMRA.Sizes(ii) < NODE_PTS_LOWER_LIMIT,
        excluded_nodes(end+1) = int64(ii-1);    % zero-based!
    else
        % IDs for vis need to be zero-based
        node_id = int64(ii-1);      % zero-based!
        node_group_name = [nodes_root int2str(node_id)];      % zero-based!

        % node_g['id'] = node['id']
        h5create(filename, [node_group_name '/id'], 1, 'Datatype', 'int64');
        h5write(filename, [node_group_name '/id'], node_id);      % zero-based!

        % if 'parent_id' in node:
        %     node_g['parent_id'] = node['parent_id']
        % cp is a list of the parent ID of each node. Root has parent ID == 0 (in Matlab's ones-based indices)
        parent_id = int64(GMRA.cp(ii)-1);      % zero-based!
        % Adopting convention (for now) that root doesn't have parent_id field
        if parent_id > -1,
            h5create(filename, [node_group_name '/parent_id'], 1, 'Datatype', 'int64');
            h5write(filename, [node_group_name '/parent_id'], parent_id);      % zero-based!
        end

        % node_g['l2Radius'] = node['l2Radius']
        h5create(filename, [node_group_name '/l2Radius'], 1, 'Datatype', 'double');
        h5write(filename, [node_group_name '/l2Radius'], GMRA.Radii(ii));

        % node_g['a'] = node['a']

        % node_g['npoints'] = node['npoints']
        h5create(filename, [node_group_name '/npoints'], 1, 'Datatype', 'int64');
        h5write(filename, [node_group_name '/npoints'], int64(GMRA.Sizes(ii)));

        d = size(GMRA.ScalFuns{ii},2);
        if d == 0,
            N = size(GMRA.ScalFuns{ii},1);
            % node_g.create_dataset('phi', data=node['phi'])
            h5create(filename, [node_group_name '/phi'], [N 1], 'Datatype', 'double');
            h5write(filename, [node_group_name '/phi'], zeros([N 1]));

            % node_g.create_dataset('sigma', data=node['sigma'])
            % NOTE: Careful, this Sigmas is normalized to 1/sqrt(n_points)...
            h5create(filename, [node_group_name '/sigma'], 1, 'Datatype', 'double');
            h5write(filename, [node_group_name '/sigma'], 0);
        else
            % node_g.create_dataset('phi', data=node['phi'])
            h5create(filename, [node_group_name '/phi'], size(GMRA.ScalFuns{ii}), 'Datatype', 'double');
            h5write(filename, [node_group_name '/phi'], GMRA.ScalFuns{ii});

            % node_g.create_dataset('sigma', data=node['sigma'])
            % NOTE: Careful, this Sigmas is normalized to 1/sqrt(n_points)...
            h5create(filename, [node_group_name '/sigma'], d, 'Datatype', 'double');
            h5write(filename, [node_group_name '/sigma'], GMRA.Sigmas{ii}(1:d));
        end
        % node_g.create_dataset('dir', data=node['dir'])
        % node_g.create_dataset('mse', data=node['mse'])

        % node_g.create_dataset('center', data=node['center'])
        h5create(filename, [node_group_name '/center'], length(GMRA.Centers{ii}), 'Datatype', 'double');
        h5write(filename, [node_group_name '/center'], GMRA.Centers{ii});

        % node_g.create_dataset('indices', data=node['indices'])
        h5create(filename, [node_group_name '/indices'], length(GMRA.PointsInNet{ii}), 'Datatype', 'int64');
        h5write(filename, [node_group_name '/indices'], int64(GMRA.PointsInNet{ii}-1));      % zero-based!

        % node_g.create_dataset('sigma2', data=node['sigma2'])
        % NOTE: Careful, this Sigmas is normalized to 1/sqrt(n_points)...
        h5create(filename, [node_group_name '/sigma2'], length(GMRA.Sigmas{ii}), 'Datatype', 'double');
        h5write(filename, [node_group_name '/sigma2'], GMRA.Sigmas{ii}.^2);
    end
end

% In high-level interface, group must exist before you can attach
% attributes to it.
h5writeatt(filename, '/full_tree', 'n_nodes', int64(n_total_nodes-length(excluded_nodes)));

%% Tree links for structure

% Second pass to put in hard links to root node and children
fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
lcpl = 'H5P_DEFAULT';
lapl = 'H5P_DEFAULT';

fprintf('Creating hard links for tree\n');
for ii = 1:length(GMRA.cp),
    % Avoid excluded nodes
    node_id = int64(ii-1);     % zero-based
    if ~any(excluded_nodes == node_id),      % zero-based!
        node_group_name = [nodes_root int2str(node_id)];      % zero-based!
        gid_node = H5G.open(fid, node_group_name);

        children_ids = (find(GMRA.cp == ii) - 1);      % zero-based!
        excluded_children = ismember(children_ids, excluded_nodes);      % zero-based!
        % Make sure there are children, and not all of them have been excluded
        if ~isempty(children_ids) && ~all(excluded_children),      % zero-based!
            % Create the new children group to put the links in
            gid_children = H5G.create(gid_node, 'children', 'H5P_DEFAULT');
            % Only loop over non-excluded children_ids
            for cc_idx = find(~excluded_children),      % zero-based!
                % Find the child node group
                child_node_name = [nodes_root int2str(int64(children_ids(cc_idx)))];      % zero-based!
                gid_child = H5G.open(fid, child_node_name);
                % Create the hard link
                H5L.create_hard(gid_child, child_node_name, gid_children, int2str(int64(children_ids(cc_idx))), lcpl, lapl);      % zero-based!
                H5G.close(gid_child);
            end
            H5G.close(gid_children);
        end
        H5G.close(gid_node);
    end
end

% Root node -- better be just one!
gid_fulltree = H5G.open(fid, '/full_tree');
root_id = (find(GMRA.cp == 0) - 1);      % zero-based!
child_node_name = [nodes_root int2str(int64(root_id))];
gid_child = H5G.open(fid, child_node_name);
% Create the hard link
H5L.create_hard(gid_child, child_node_name, gid_fulltree, 'tree_root', lcpl, lapl);
H5G.close(gid_child);
H5G.close(gid_fulltree);

H5F.close(fid);

end
