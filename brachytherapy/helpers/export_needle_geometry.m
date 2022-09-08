% export needle entry point and target point for the deformation simulation 

clear all 
close all

%% load seed data
patient = 4;
geometry_path = ['~/Results/2022_07_06 needle geometries/Pat' num2str(patient), '/'];
dose_path = ['~/Results/2022_07_07 dose planning/Pat' num2str(patient), '/2.SA/p400,r30,u1/'];
load([geometry_path, 'tplan_full.mat'], "full_tplan")
load([dose_path, 'result.mat'], "resultGUI")
load([geometry_path, 'moved_needles.mat'], "moved_needles")


%% extract points
weights = resultGUI.w;
used_points = [];
for i = 1:numel(weights)
    if(weights(i) == 1)
        used_points = [used_points; full_tplan{i,2}'];
    end
end

%% extract seeds for each needle
n = 1;
ind = [];
for j = 1:(size(used_points,1)-1)
    current_z = used_points(j,3);
    if (used_points(j+1,3)>used_points(j,3) && ...
            ((used_points(j+1,1)-used_points(j,1))^2 + (used_points(j,2)-used_points(j+1,2))^2) < 5)
        ind = [ind, j];
    else
        ind = [ind,j];
        % check if points are in simular x-y position as cog
        needle_found = false;
        while (~needle_found)     
            if (((used_points(ind(1),1)-moved_needles{n,1}.center_of_gravity(1))^2 + ...
                (used_points(ind(1),2)-moved_needles{n,1}.center_of_gravity(2))^2) < 5)
                moved_needles{n,1}.used_points = used_points(ind,:);
                n = n+1;
                ind = [];
                needle_found = true;
            else
                n = n+1;
            end
        end
    end
end

%% get target point and angle and delete needles which are not used
% equation for angle check Labbook, 13.07.22
needles_to_delete = [];
b = [0,1,0];
for i = 1:numel(moved_needles)
    if (isfield(moved_needles{i,1},'used_points'))
        moved_needles{i,1}.target_point = moved_needles{i,1}.used_points(end,:);
        moved_needles{i,1}.target_point_ortho = [moved_needles{i,1}.template_entry_point(1), ...
            moved_needles{i,1}.template_entry_point(2),moved_needles{i,1}.target_point(3)];
        k = moved_needles{i,1}.target_point - moved_needles{i,1}.target_point_ortho;
        if (norm(k) == 0)
            moved_needles{i}.alpha = 0;
        else
            moved_needles{i}.alpha = acosd((k*b')/(norm(k)*norm(b)));
        end
    else
        needles_to_delete = [needles_to_delete, i];
    end
end
moved_needles(needles_to_delete) = [];

%% plot for control
figure
hold on
for i = 1:numel(moved_needles)
    plot3(moved_needles{i}.template_entry_point(1),moved_needles{i}.template_entry_point(:,2),moved_needles{i}.template_entry_point(:,3), '-o', 'Color', [1,0,1],'MarkerFaceColor',[1, 1-i/17, i/17],'MarkerEdgeColor',[1, 1-i/17, i/17]);
    quiver3(moved_needles{i}.template_entry_point(1),moved_needles{i}.template_entry_point(:,2),moved_needles{i}.template_entry_point(:,3), ...
        moved_needles{i}.needle_direction(1),moved_needles{i}.needle_direction(:,2),moved_needles{i}.needle_direction(:,3), 100)
     if (isfield(moved_needles{i,1},'used_points'))
         plot3(moved_needles{i}.target_point(1),moved_needles{i}.target_point(:,2),moved_needles{i}.target_point(:,3), '-o', 'Color', [1,0,1],'MarkerFaceColor',[1, 1-i/17, i/17],'MarkerEdgeColor',[1, 1-i/17, i/17]);
     end
end

%% export
JSONFILE_name = sprintf('configNeedles.json'); 
fid = fopen(JSONFILE_name,'w'); 
Tip_force_magnitude = 140;
Width = 3e-3;
Beam_length = 12.5e-3;
Number_of_beam_elements = 16;
s = struct();
for i = 1:size(moved_needles,1)
    s(i).Instrument_ID = i;
    s(i).Entry_point = 1e-3*( moved_needles{i}.template_entry_point + [50,70,85]);  % coordinate transform to FEM grid
    s(i).Target_point = 1e-3*( moved_needles{i}.target_point_ortho + [50,70,85]);   % coordinate transform to FEM grid        
    s(i).Rotation = moved_needles{i}.alpha; 
    s(i).Tip_force_magnitude = Tip_force_magnitude; 
    s(i).Width = Width;
    s(i).Beam_length = Beam_length;
    s(i).Number_of_beam_elements = Number_of_beam_elements;
end
n = struct("Needles", s);
encodedJSON = jsonencode(n,PrettyPrint=true);
fprintf(fid, encodedJSON); 
 
fclose('all'); 