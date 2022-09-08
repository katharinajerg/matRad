% reduce grid size in vtu unstructured grid
deformationDataPath = '/home/kjerg/Results/2022_08_18 tissue elasticity/Pat1/46800_140_results_physical_domain_10.vtu';
deformationDataPathOut = '/home/kjerg/Results/2022_08_18 tissue elasticity/Pat1/46800_140_results_physical_domain_1_reduced.vtu';


%% import data
tissue(1) = struct(); 
imported_tissue = vtkRead(deformationDataPath);
tissue.verticies = 1000*imported_tissue.points;
tissue.deformation_field = 1000*imported_tissue.pointData.displacement;
clear imported_tissue

%%
shifts = [-50,-70,-85];
tissue.verticies(:,1) = tissue.verticies(:,1) + shifts(1); 
tissue.verticies(:,2) = tissue.verticies(:,2) + shifts(2);
tissue.verticies(:,3) = tissue.verticies(:,3) + shifts(3);

%% Remove verticies outside of physical domain (phase field < 0.5). -- does not work, because phase field is not saved
%ind_outside_physical_domain = find(all(tissue.phase_field < 0.5,2));
%tissue.verticies(ind_outside_physical_domain,:) = [];
%tissue.deformation_field(ind_outside_physical_domain,:) = [];
%tissue.phase_field(ind_outside_physical_domain) = [];

%% instead: remove all deformations = zero, besides back wall, which is zero due to boundary conditions
ind_zero_deformation = find(all(tissue.deformation_field < 1e-6,2));
ind_back_wall = find(tissue.verticies(:,3) > ((150+shifts(3))-1e-6));
ind_no_phase_field = ismember(ind_zero_deformation, ind_back_wall);
ind_phase_field = ind_zero_deformation;
ind_phase_field(ind_no_phase_field) = [];
tissue.verticies(ind_phase_field,:) = [];
tissue.deformation_field(ind_phase_field,:) = [];

%% Remove doubled occurance of verticies.
[less_vert, ia, ic] = unique(tissue.verticies, 'stable', 'rows');
tissue.verticies = tissue.verticies(ia,:);
tissue.deformation_field = tissue.deformation_field(ia,:);

%% delete points

%% plot
close all
tissue_plot = tissue;
ind = find(tissue.verticies(:,2) > (-31) | tissue.verticies(:,2) < (-32));
tissue_plot.verticies(ind,:)=[];
tissue_plot.deformation_field(ind,:)=[];
ind = floor(linspace(1,size(tissue_plot.verticies,1),round(0.8*size(tissue_plot.verticies,1)))');
tissue_plot.verticies(ind,:)=[];
tissue_plot.deformation_field(ind,:)=[];
figure
% scatter3(tissue_plot.verticies(:,1), tissue_plot.verticies(:,2), tissue_plot.verticies(:,3),'.')
% hold on 
quiver3(tissue_plot.verticies(:,1), tissue_plot.verticies(:,2), tissue_plot.verticies(:,3),...
    tissue_plot.deformation_field(:,1), tissue_plot.deformation_field(:,2), tissue_plot.deformation_field(:,3), 'k');
figure
scatter3(tissue_plot.verticies(:,1), tissue_plot.verticies(:,2), ...
    tissue_plot.verticies(:,3),'.')


