% ========================================================================
% GMSH_IMPORT: Import a .m-file generated by GMSH (containing a 'msh' struct)
%              and convert it into a structured format for body-wise processing.
%
% Usage:
%   satellite = gmsh_import('C:/path/to/mesh/full_satellite.m');
%
% The input file must define a variable 'msh' with fields:
%   - msh.POS:     Node positions (n_nodes × 3)
%   - msh.TRIANGLES: Element definitions (n_elements × 4), where
%                    columns 1:3 = node indices of each triangle
%                    column 4   = body ID
%
% Output:
%   body_data struct with:
%     - vertices_CAD:  cell array of size {n_bodies×1}, each cell is 3×3×n_triangles
%     - centroids_CAD: cell array {n_bodies×1}, each 3×n matrix of triangle centroids
%     - normals_CAD:   cell array {n_bodies×1}, each 3×n matrix of normal vectors
%     - areas_CAD:     cell array {n_bodies×1}, each 1×n vector of triangle areas
% ========================================================================

function body_data = import_gmsh(relative_path)

    % Save current working directory to restore it later
    old_dir = pwd;

    % Separate file path, name and extension
    [filepath, filename, ext] = fileparts(relative_path);

    % If a folder is specified, switch into it
    if ~isempty(filepath)
        cd(filepath);
    end

    clear msh;  % Ensure no leftover msh variable from workspace

    full_path = fullfile(filepath, filename+ ext);
    % Check if the file exists
    if ~exist(full_path, 'file')
        error('File not found: %s', full_path);
    end

    % Run the GMSH-exported .m file (should define 'msh')
    run(full_path);


    % Validate msh struct existence
    if ~exist('msh', 'var')
        error('The GMSH-exported file does not define a "msh" structure.');
    end

    % Get all unique body IDs from the 4th column of the triangle data
    body_ids = unique(msh.TRIANGLES(:,4));

    % Extract vertex coordinates for each body
    % For each body:
    % - collect all triangle node indices (3 per triangle)
    % - use these to extract XYZ coordinates from msh.POS
    % - reshape to a 3×3×n array (3 nodes × 3 coordinates × n triangles)
    vertices_B = arrayfun(@(b) ...
        reshape(msh.POS(reshape(msh.TRIANGLES(msh.TRIANGLES(:,4) == b, 1:3)', [], 1), :)', ...
                [3, 3, nnz(msh.TRIANGLES(:,4) == b)]), ...
        body_ids, 'UniformOutput', false);

    % Return to original working directory
    cd(old_dir);

    % Create output structure
    body_data = struct('vertices_CAD', vertices_B);

end

