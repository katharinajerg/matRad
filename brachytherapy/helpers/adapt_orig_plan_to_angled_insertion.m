clear all
close all

% Pat4
path = '~/Daten/';
D = '~/thindrives/Brachy18_02/dicom-dict-iotp.txt';
info_pl = dicominfo([path,'PL001.dcm'],'dictionary',D);

% source locations needle insertion
seeds = struct2array(info_pl.ApplicationSetupSequence);
tplan = cell(length(seeds),3);

for i = 1:length(seeds)
    tplan{i,1} = seeds(i).ApplicationSetupNumber;
    tplan{i,2} = seeds(i).ChannelSequence.Item_1.BrachyControlPointSequence.Item_1.ControlPoint3DPosition;
    tplan{i,3} = seeds(i).ChannelSequence.Item_1.BrachyControlPointSequence.Item_1.ControlPointOrientation;
end

% MR Dicom
info_mr = dicominfo([path,'MR001.dcm']);
mr = dicomread(info_mr);
grid = zeros(size(mr));
grid_points = zeros(size(mr));
grid(mr == 229) = 1;
grid_coordinates = [0,0];

% filter middle of grid marker
for i = 2: size(grid_points,1)-1
    for j = 2:size(grid_points,2)-1
        if (grid(i,j) == 1 && grid(i-1,j) == 1 && grid(i+1,j) == 1 && grid(i,j-1) == 1 && grid(i,j+1) == 1 ...
                && grid(i-1,j-1) == 0 && grid(i-1,j+1) == 0 && grid(i+1,j-1) == 0 && grid(i+1,j+1) == 0)
            grid_points(i,j) = 1;
            grid_coordinates = [grid_coordinates; [j,i]];
        end
    end
end
grid_coordinates(1,:) = [];

% transform grid coordinates from voxels to patient coordinates 
pixel_spacing = info_mr.PixelSpacing;
image_position = info_mr.ImagePositionPatient;
image_orientation = info_mr.ImageOrientationPatient;
grid_coordinates_patient = zeros(size(grid_coordinates));

% equation from: https://dicom.nema.org/medical/Dicom/2016b/output/chtml/part03/sect_C.7.6.2.html#sect_C.7.6.2.1.1
% grid coordinates -1, because in equation first column/row index is zero.
for i = 1:size(grid_coordinates,1)
    grid_coordinates_patient(i,1) = image_orientation(1)*pixel_spacing(1)*(grid_coordinates(i,1)-1) + ...
                                    image_orientation(4)*pixel_spacing(2)*(grid_coordinates(i,2)-1) + image_position(1);
    grid_coordinates_patient(i,2) = image_orientation(2)*pixel_spacing(1)*(grid_coordinates(i,1)-1) + ...
                                    image_orientation(5)*pixel_spacing(2)*(grid_coordinates(i,2)-1) + image_position(2);
end


% extract needle paths 
needles = struct('points', {}, 'number_of_seeds', {});
point_array(1,:) = tplan{1,2};
needle_index = 1;
point_index = 2;
for i = 2:length(tplan)
    diff = norm(tplan{i,2}(1:2)-tplan{i-1,2}(1:2));
    if (abs(diff) < 2 )
        point_array(point_index,:) = tplan{i,2};
        point_index = point_index + 1;
    else
        needles{needle_index,1}.points = point_array;
        needles{needle_index,1}.number_of_seeds = point_index-1;
        point_array = [];
        point_array(1,:) = tplan{i,2};
        point_index = 2;
        needle_index = needle_index + 1;
    end
end
needles{needle_index,1}.points = point_array;
needles{needle_index,1}.number_of_seeds = point_index-1;

% if seeds have been added later, they are recognized as seperate needle
% therefore concatinate needles with corresponding x-y-values of seed
% position
js = [];
for i = 1:numel(needles)
    for j = i+1:numel(needles)
        pos_i = needles{i}.points(1,1:2);
        pos_j = needles{j}.points(1,1:2);
        diff = norm(pos_j - pos_i);
        if diff < 2
            needles{i}.points(end+1,:) = needles{j}.points(:,:);
            needles{i}.number_of_seeds = needles{i}.number_of_seeds + needles{j}.number_of_seeds;
            [~, ind] = sort(needles{i}.points(:,3),1,'descend');
            needles{i}.points = needles{i}.points(ind, :);
            js = [js, j];
        end
    end
end
needles(js) = [];

% find suitable entry point and calculate center of gravity of seeds
grid_coordinates_patient(:,3) = -85;
figure 
hold on 
scatter3(grid_coordinates_patient(:,1), grid_coordinates_patient(:,2), grid_coordinates_patient(:,3), 'r+');

all_ind = [];
for i = 1:length(needles)
    min_dist = 100;
    ind = 0;
    for j = 1:size(grid_coordinates_patient,1)
        dist = norm(grid_coordinates_patient(j,1:2)-needles{i}.points(1,1:2));
        if (dist < min_dist && ~sum(all_ind == j))
            min_dist = dist;
            ind = j;
        end
    end
    plot3(needles{i}.points(:,1),needles{i}.points(:,2),needles{i}.points(:,3), '-o', 'Color', [0,0,1],'MarkerFaceColor',[0, 1-i/needle_index, i/needle_index],'MarkerEdgeColor',[0, 1-i/needle_index, i/needle_index]);
    needles{i}.template_entry_point = grid_coordinates_patient(ind,:);
    needles{i}.center_of_gravity = [sum(needles{i}.points(:,1))/needles{i}.number_of_seeds, ...
                                    sum(needles{i}.points(:,2))/needles{i}.number_of_seeds, ...
                                    sum(needles{i}.points(:,3))/needles{i}.number_of_seeds];
    needles{i}.needle_direction = [0,0,1];
    all_ind = [all_ind, ind];
end

% for all possible seed points on needle, fill in remaining points
full_needles = needles;
for i = 1:numel(full_needles)
    min_z_value = full_needles{i}.points(end,3);
    max_z_value = full_needles{i}.points(1,3);
    % check if min is acutally smaller or equal to max
    if min_z_value > max_z_value
        temp = min_z_value;
        min_z_value = max_z_value;
        max_z_value = temp;
    end
    % 5 mm margin in each direction 
    min_z_value = min_z_value - 10;
    max_z_value = max_z_value + 10;
    z_values = min_z_value:5:max_z_value;
    points = repmat(full_needles{i}.points(1,1:2), [numel(z_values), 1]);
    points = [points, z_values'];
    full_needles{i}.points = points;
    full_needles{i}.number_of_seeds = size(points, 1);
end
full_needles_orig = full_needles;


% calculate new needle direction and corrected seed positions for both
% needle configurations
moved_needles = needles;
for i = 1:numel(needles)
    direction = needles{i}.center_of_gravity - needles{i}.template_entry_point;
    moved_needles{i}.needle_direction = direction/norm(direction);
    for s = 1:moved_needles{i}.number_of_seeds
        dist_to_cog = needles{i}.points(s,:) - needles{i}.center_of_gravity;
        if dist_to_cog(3) > 0
            moved_needles{i}.points(s,:) = needles{i}.center_of_gravity + norm(dist_to_cog) * moved_needles{i}.needle_direction;
        else
            moved_needles{i}.points(s,:) = needles{i}.center_of_gravity - norm(dist_to_cog) * moved_needles{i}.needle_direction;
        end
    end
    moved_needles{i}.points(end+1,:) = needles{i}.template_entry_point;
    plot3(moved_needles{i}.points(:,1),moved_needles{i}.points(:,2),moved_needles{i}.points(:,3), '-o', 'Color', [1,0,1],'MarkerFaceColor',[1, 1-i/needle_index, i/needle_index],'MarkerEdgeColor',[1, 1-i/needle_index, i/needle_index]);

    full_needles{i}.needle_direction = direction/norm(direction);
    for s = 1:full_needles{i}.number_of_seeds
        dist_to_cog = full_needles{i}.points(s,:) - needles{i}.center_of_gravity;
        if dist_to_cog(3) > 0
            full_needles{i}.points(s,:) = needles{i}.center_of_gravity + norm(dist_to_cog) * full_needles{i}.needle_direction;
        else
            full_needles{i}.points(s,:) = needles{i}.center_of_gravity - norm(dist_to_cog) * full_needles{i}.needle_direction;
        end
    end
    full_needles{i}.points(end+1,:) = needles{i}.template_entry_point;
end

hold off
title('needle paths (needle insertion)')

% write data to tplan and save 
moved_tplan = tplan;
ind = 0;
for n = 1 : numel(needles)
    for s = 1:moved_needles{n}.number_of_seeds
        ind = ind + 1;
        moved_tplan{ind,2}(1,1) = moved_needles{n}.points(s,1);
        moved_tplan{ind,2}(2,1) = moved_needles{n}.points(s,2);
        moved_tplan{ind,2}(3,1) = moved_needles{n}.points(s,3);
        moved_tplan{ind,3}(1,1) = moved_needles{n}.needle_direction(1);
        moved_tplan{ind,3}(2,1) = moved_needles{n}.needle_direction(2);
        moved_tplan{ind,3}(3,1) = moved_needles{n}.needle_direction(3);
    end
end
full_tplan = cell(1,3);
ind = 0;
for n = 1 : numel(needles)
    for s = 1:full_needles{n}.number_of_seeds
        ind = ind + 1;
        full_tplan{ind,1}(1) = ind;
        full_tplan{ind,2}(1,1) = full_needles{n}.points(s,1);
        full_tplan{ind,2}(2,1) = full_needles{n}.points(s,2);
        full_tplan{ind,2}(3,1) = full_needles{n}.points(s,3);
        full_tplan{ind,3}(1,1) = full_needles{n}.needle_direction(1);
        full_tplan{ind,3}(2,1) = full_needles{n}.needle_direction(2);
        full_tplan{ind,3}(3,1) = full_needles{n}.needle_direction(3);
    end
end
full_tplan_orig = cell(1,3);
ind = 0;
for n = 1 : numel(needles)
    for s = 1:full_needles{n}.number_of_seeds
        ind = ind + 1;
        full_tplan_orig{ind,1}(1) = ind;
        full_tplan_orig{ind,2}(1,1) = full_needles_orig{n}.points(s,1);
        full_tplan_orig{ind,2}(2,1) = full_needles_orig{n}.points(s,2);
        full_tplan_orig{ind,2}(3,1) = full_needles_orig{n}.points(s,3);
        full_tplan_orig{ind,3}(1,1) = full_needles_orig{n}.needle_direction(1);
        full_tplan_orig{ind,3}(2,1) = full_needles_orig{n}.needle_direction(2);
        full_tplan_orig{ind,3}(3,1) = full_needles_orig{n}.needle_direction(3);
    end
end
save("tplan.mat", "moved_tplan")
save("tplan_full.mat", "full_tplan")
save("tplan_full_orig.mat", "full_tplan_orig")
save("tplan_orig.mat", "tplan")

% determine bool vector with 1 for each seed in full_tplan with is occupied
% by a seed in t_plan
moved_pos = reshape(cell2mat(moved_tplan(:,2)),[3,size(moved_tplan,1)])';
full_pos = reshape(cell2mat(full_tplan(:,2)),[3,size(full_tplan,1)])';
supp = ismember(full_pos,moved_pos,'rows');
save("supp.mat", "supp")

