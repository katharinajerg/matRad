% Katharina Jerg, June 2023
% This script is used to evaluate the results of the differentiation of the
% dose contraints regarding the seed positions.

%% prepare data and calculate distances to surfaces
clear all
close all

% patient data sets
patientIds = 1:35;
targetNames = {'P_V100', 'P_D90', 'R_D2cc', 'U_D30'};
data = cell(size(patientIds,2),1);  

% load data
for i = 1:size(patientIds,2)
    id = patientIds(i)+1000;
    path = '..\BRACHYTHERAPY_data\evaluation\';
    gradval_U_D30 = load([path, 'U_D30\gradval_',num2str(id),'.mat']).gradval;
    fval_U_D30 = load([path, 'U_D30\fval_',num2str(id),'.mat']).fval;
    gradval_P_D90 = load([path, 'P_D90\gradval_',num2str(id),'.mat']).gradval;
    fval_P_D90 = load([path, 'P_D90\fval_',num2str(id),'.mat']).fval;
    gradval_P_V100 = load([path, 'P_V100\gradval_',num2str(id),'.mat']).gradval;
    fval_P_V100 = load([path, 'P_V100\fval_',num2str(id),'.mat']).fval;
    gradval_R_D2cc = load([path, 'R_D2cc\gradval_',num2str(id),'.mat']).gradval;
    fval_R_D2cc = load([path, 'R_D2cc\fval_',num2str(id),'.mat']).fval;

    % fill data struct with information 
    for t = 1:length(targetNames)
        targetName = targetNames{t};
        gradval_name  = [ 'gradval_',targetName];
        val_name  = [ 'fval_',targetName];
        fieldname_grad  = [targetName, '_gradients'];
        fieldname_gradMag = [targetName, '_gradMagnitudes'];
        fieldname_val = [targetName];
        fieldname_mean = [targetName, '_mean_gradMagnitudes'];
        fieldname_max = [targetName, '_max_gradMagnitudes'];
        gradPerSeed = reshape(extractdata(eval(gradval_name)), [3,size(eval(gradval_name),2)/3]);
        magnitudePerSeed = vecnorm(gradPerSeed);
        data{i}.(fieldname_grad) = gradPerSeed;
        data{i}.(fieldname_gradMag) = magnitudePerSeed;
        data{i}.(fieldname_val) = extractdata(eval(val_name));
        data{i}.(fieldname_mean) = median(magnitudePerSeed);
        data{i}.(fieldname_max) = max(magnitudePerSeed);
    end
        data{i}.id = id;

    % prostate information
    path = ['..\BRACHYTHERAPY_data\',num2str(id-1000),'\IntraOp\IntraOp\'];
    pathStructureSet = [path, 'SS001.dcm'];
    pathImg = [path, 'MR001.dcm'];
    [cst, ct] = matRad_importDicomUSStructureSet(pathStructureSet,pathImg);
    cube = zeros(size(ct.cube{1}));
    cube(cst{1,4}{1}) = 1;
    data{i}.prostateVolume = sum(cube(:));

    %   extend in x, y, z
    line_x = sum(cube,3);
    line_x = sum(line_x,1); % x is dim = 2
    x_min =  ct.x(find(line_x, 1, 'first'));
    x_max =  ct.x(find(line_x, 1, 'last'));
    data{i}.xExtendProstate = [x_min, x_max];  

    line_y = sum(cube,3);
    line_y = sum(line_y,2); % y is dim = 1
    y_min =  ct.y(find(line_y, 1, 'first'));
    y_max =  ct.y(find(line_y, 1, 'last'));
    data{i}.yExtendProstate = [y_min, y_max];

    line_z = sum(cube,1);
    line_z = sum(line_z,2);
    z_min =  ct.z(find(line_z, 1, 'first'));
    z_max =  ct.z(find(line_z, 1, 'last'));
    data{i}.zExtendProstate = [z_min, z_max];

    %   center of mass
    [X, Y, Z] = meshgrid(ct.x, ct.y, ct.z);
    x_com = sum(sum(sum(X.*cube))) / sum(cube(:));
    y_com = sum(sum(sum(Y.*cube))) / sum(cube(:));
    z_com = sum(sum(sum(Z.*cube))) / sum(cube(:));
    data{i}.centerOfMass = [x_com, y_com, z_com];

    % seed positions
    load(['..\BRACHYTHERAPY_data\',num2str(id-1000),'\IntraOp\IntraOp\tplan_orig.mat']);

    % for each seed get distance to urethra, rectum and prostate surface
    factor = 1;
    urethra_cube = zeros(size(ct.cube{1}));
    urethra_cube(cst{2,4}{1}) = 1;
    urethra_cube = imresize3(urethra_cube, factor);
    surface_urethra = isosurface(X, Y, Z,urethra_cube, 0.5);
    [distance_map_urethra, idx_urethra] = bwdist(urethra_cube);

    rectum_cube = zeros(size(ct.cube{1}));
    rectum_cube(cst{3,4}{1}) = 1;
    rectum_cube= imresize3(rectum_cube, factor);
    surface_rectum = isosurface(X, Y, Z,rectum_cube, 0.5);
    [distance_map_rectum, idx_rectum] = bwdist(rectum_cube);

    prostate_cube = double(~cube);
    prostate_cube= imresize3(prostate_cube, factor);

    % close prostate surface
    X = cat(3,X(:,:,1),X);
    X = cat(3,X,X(:,:,end));
    Y = cat(3,Y(:,:,1),Y);
    Y = cat(3,Y,Y(:,:,end));
    z_1_val = Z(1,1,1) - (Z(1,1,2)-Z(1,1,1));
    z_2_val = Z(1,1,end) + (Z(1,1,2)-Z(1,1,1));
    Z = cat(3,ones(size(Z,1,2))*z_1_val, Z);
    Z = cat(3,Z,ones(size(Z,1,2))*z_2_val);
    prostate_cube = cat(3,ones(size(prostate_cube,1,2)), prostate_cube);
    prostate_cube = cat(3,prostate_cube, ones(size(prostate_cube,1,2)));

    surface_prostate = isosurface(X, Y, Z, prostate_cube, 0.5);
    [distance_map_prostate, idx_prostate] = bwdist(prostate_cube);

    voxel_size = factor*[ct.resolution.x, ct.resolution.y, ct.resolution.z];

    % convert to column vector
    cube_size = size(prostate_cube);
    ind_cube = zeros(size(prostate_cube));
    for x_i = 1:cube_size(1)
        for y_i = 1:cube_size(2)
            for z_i = 1:cube_size(3)
                ind_cube(x_i,y_i,z_i) = cube_size(1)*cube_size(2)*(z_i-1) + cube_size(1)*(y_i-1) + x_i;
            end
        end
    end

    data{i}.n_prostate = [];
    data{i}.n_rectum = [];
    data{i}.n_urethra = [];

    for s = 1:size(tplan,1)
    
        point_coordinates = cell2mat(tplan(s,2))';

        % get distances to and points on triangular surface
        % prostate
        [ distance_point_prostate, point_on_prostate, ~, ~, ~, ~] = point2trimesh('Faces', surface_prostate.faces, 'Vertices', surface_prostate.vertices,  'QueryPoints', point_coordinates );
        n_prostate = (point_coordinates - point_on_prostate) / norm((point_coordinates - point_on_prostate));
        data{i}.n_prostate = [data{i}.n_prostate, n_prostate'];
        tplan(s,4) = mat2cell(abs(distance_point_prostate),1); % column 4: prostate
        % rectum
        [ distance_point_rectum, point_on_rectum, ~, ~, ~, ~] = point2trimesh('Faces', surface_rectum.faces, 'Vertices', surface_rectum.vertices,  'QueryPoints', point_coordinates );
        n_rectum= (point_coordinates - point_on_rectum) / norm((point_coordinates - point_on_rectum));
        data{i}.n_rectum= [data{i}.n_rectum, n_rectum'];
        tplan(s,5) = mat2cell(abs(distance_point_rectum),1); % column 5: rectum
        % urethra
        [ distance_point_urethra, point_on_urethra, ~, ~, ~, ~] = point2trimesh('Faces', surface_urethra.faces, 'Vertices', surface_urethra.vertices,  'QueryPoints', point_coordinates );
        n_urethra= (point_coordinates - point_on_urethra) / norm((point_coordinates - point_on_urethra));
        data{i}.n_urethra = [data{i}.n_urethra, n_urethra'];
        tplan(s,6) = mat2cell(abs(distance_point_urethra),1); % column 6: rectum

    end

    data{i}.tplan = tplan;

    % seed activity
    info_pl = dicominfo(['..\BRACHYTHERAPY_data\',num2str(i),'\IntraOp\IntraOp\PL001.dcm']);
    data{i}.activity = info_pl.SourceSequence.Item_1.ReferenceAirKermaRate;
end

clearvars -except data patientIds targetNames

% evaluate data
% extract sensitive seeds
evaluation_data = cell(1,1);
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_gradMag = [targetName, '_gradMagnitudes'];
    evaluation_data{1}.(fieldname_gradMag) = [];
end

for i=1:size(patientIds,2)
    for t = 1:length(targetNames)
        targetName = targetNames{t};
        fieldname_gradMag = [targetName, '_gradMagnitudes'];
        evaluation_data{1}.(fieldname_gradMag) = [evaluation_data{1}.(fieldname_gradMag); data{i}.(fieldname_gradMag)'];
    end
end
% Seed considered sensitive are the seeds which have a gradient magnitude
% which is within the largest 5% of the entire data.
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_gradMagnitudes = [targetName, '_gradMagnitudes'];
    fieldname_gradMagnitudesSorted = [targetName, '_gradMagnitudesSorted'];
    fieldname_threshold = [targetName, '_threshold'];
    fieldname_numOutliers = [targetName, '_numOutliersPerPatient'];
    evaluation_data{1}.(fieldname_numOutliers) = [];
    [evaluation_data{1}.(fieldname_gradMagnitudesSorted),ind] = sort(evaluation_data{1}.(fieldname_gradMagnitudes));
    evaluation_data{1}.(fieldname_threshold) = evaluation_data{1}.(fieldname_gradMagnitudesSorted)(floor(size(evaluation_data{1}.(fieldname_gradMagnitudesSorted),1)*0.95));
end

for i=1:size(patientIds,2)
    numOutliers = struct();
    for t = 1:length(targetNames)
        targetName = targetNames{t};
        numOutliers.(targetName) = 0;
        outlierIdx.(targetName) = [];
    end
    for s = 1:size(data{i}.tplan,1)
        seedPosition = data{i}.tplan{s,2};
        % if(ismember(s, data{i}.outlierIdxs)) 
        boolOutliers = 0;
        for t = 1:length(targetNames)
            targetName = targetNames{t};
            fieldname_gradMagnitudes = [targetName, '_gradMagnitudes'];
            fieldname_threshold = [targetName, '_threshold'];
            fieldname_outliers = [targetName, '_outlierIdx'];

            if(data{i}.(fieldname_gradMagnitudes)(s) > evaluation_data{1}.(fieldname_threshold))
                numOutliers.(targetName) = numOutliers.(targetName) + 1;
                outlierIdx.(targetName) = [outlierIdx.(targetName), s];
            end
            data{i}.(fieldname_outliers) = outlierIdx.(targetName);
        end
    end
    for t = 1:length(targetNames)
        targetName = targetNames{t};
        fieldname_numOutlierspP = [targetName, '_numOutliersPerPatient'];
        evaluation_data{1}.(fieldname_numOutlierspP) = [evaluation_data{1}.(fieldname_numOutlierspP), numOutliers.(targetName)];    
    end
end


%% plot 3D views
% load data
for i = 29
    % for t = 1:length(targetNames)
    %     targetName = targetNames{t};
    %     fieldname_gradMag = [targetName, '_gradMagnitudes'];
       
        id = patientIds(i)+1000;
        % prostate information
        path = ['..\BRACHYTHERAPY_data\',num2str(id-1000),'\IntraOp\IntraOp\'];
        pathStructureSet = [path, 'SS001.dcm'];
        pathImg = [path, 'MR001.dcm'];
        [cst, ct] = matRad_importDicomUSStructureSet(pathStructureSet,pathImg);
        [X, Y, Z] = meshgrid(ct.x, ct.y, ct.z);
        cube = zeros(size(ct.cube{1}));
        cube(cst{1,4}{1}) = 1;
        data{i}.prostateVolume = sum(cube(:));
    
        % seed positions
        load(['..\BRACHYTHERAPY_data\',num2str(id-1000),'\IntraOp\IntraOp\tplan_orig.mat']);

        % prepare figure
        width = 800;
        height = 800;
        fontsize = 28;
        fontname = 'Times New Roman';
        figure('Color', 'w', 'Position', [100, 100, width, height]);
        hold on 
        view(3)
        ax = gca; 
        xlim([-27 27])
        ylim([-60 0])
        zlim([-40 0])
        axis off

        % plot slice (dummy, which is out of visible range in order to get
        % axis settings from wrap)
        camera_pos = ax.CameraPosition;
        camera_target = ax.CameraTarget;
        slice_num = 8;
        warp(X(:,:,slice_num),Y(:,:,slice_num),20*ones(size(X(:,:,slice_num))),zeros(size(X(:,:,slice_num))));
        ax.CameraPosition = camera_pos;
        ax.CameraTarget = camera_target;
        ax.View = [-37.5,30];

        % plot seeds halos
        ax2 = axes;
        warp(X(:,:,slice_num),Y(:,:,slice_num),-50*ones(size(X(:,:,slice_num))),zeros(size(X(:,:,slice_num))));
        ax2.Color = 'none';
        ax2.View = ax.View;
        ax2.XLim = ax.XLim;
        ax2.YLim = ax.YLim;
        ax2.ZLim = ax.ZLim;
        ax2.CameraPosition = camera_pos;
        ax2.CameraTarget = camera_target;
        axis off
        %S = 80;
        hold on 
        %targetNamesHalo = {'P_D90', 'R_D2cc', 'U_D30'};
        targetNamesHalo = {'R_D2cc'};
        gradients = zeros(size(tplan,1), length(targetNamesHalo));
        for s = 1:size(tplan,1)
            for t = 1:length(targetNamesHalo)
                targetName = targetNamesHalo{t};
                fieldname_gradMag = [targetName, '_gradMagnitudes'];
                point_coordinates = cell2mat(tplan(s,2))';
                gradients(s,t) = data{i}.(fieldname_gradMag)(s);
            end
        end
        allGradients = sum(gradients,2);
        allInverseGradients = 1.5*ones(size(allGradients))./allGradients;
        %allInverseGradients = 3*ones(size(allGradients))./allGradients;
        for s = 1:size(tplan,1)
            point_coordinates = cell2mat(tplan(s,2))';
            [Xcirc,Ycirc,Zcirc] = ellipsoid(point_coordinates(1),point_coordinates(2),point_coordinates(3),allInverseGradients(s),allInverseGradients(s),0.001);
            x_col = Xcirc-point_coordinates(1);
            y_col = Ycirc-point_coordinates(2);
            r = sqrt(x_col.^2 + y_col.^2)/allInverseGradients(s);
            r_seed = 0.5;
            [Xseed,Yseed,Zseed] = ellipsoid(point_coordinates(1),point_coordinates(2),point_coordinates(3),r_seed, r_seed, 2);
            if allInverseGradients(s) < 5
                surf(Xcirc,Ycirc,Zcirc,r,'FaceAlpha',0.5,'EdgeColor','none');
                surf(Xseed,Yseed,Zseed,'FaceColor',[0 0 0], 'FaceAlpha',0.5,'EdgeColor','none');
            else
                surf(Xseed,Yseed,Zseed,'FaceColor',[0 0.9 0.1], 'FaceAlpha',0.5,'EdgeColor','none');
            end
            %scatter3(point_coordinates(1),point_coordinates(2),point_coordinates(3),S, data{i}.(fieldname_gradMag)(s), 'filled')
        end      
        cmap = flipud(hsv);
        cmap(1:156, :) = repmat(cmap(156,:), 156,1);
        colormap(cmap);
        m = length(cmap);
        %cindex = fix((data{i}.(fieldname_gradMag)-cmin)/(cmax-cmin)*m)+1;

        % plot organ surfaces
        factor = 1;
        urethra_cube = zeros(size(ct.cube{1}));
        urethra_cube(cst{2,4}{1}) = 1;
        urethra_cube = imresize3(urethra_cube, factor);
        surface_urethra = isosurface(X, Y, Z,urethra_cube, 0.5);
        p = patch(surface_urethra);
        set(p,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor','none');  
        camlight;
        lighting gouraud;
    
        rectum_cube = zeros(size(ct.cube{1}));
        rectum_cube(cst{3,4}{1}) = 1;
        rectum_cube= imresize3(rectum_cube, factor);
        surface_rectum = isosurface(X, Y, Z,rectum_cube, 0.5);
    
        p = patch(surface_rectum);
        set(p,'FaceColor',[0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor','none');  
        set(p,'EdgeColor','none');
        camlight;
        lighting gouraud;
    
        prostate_cube = double(~cube);
        prostate_cube= imresize3(prostate_cube, factor);

        % close prostate surface
        X = cat(3,X(:,:,1),X);
        X = cat(3,X,X(:,:,end));
        Y = cat(3,Y(:,:,1),Y);
        Y = cat(3,Y,Y(:,:,end));
        z_1_val = Z(1,1,1) - (Z(1,1,2)-Z(1,1,1));
        z_2_val = Z(1,1,end) + (Z(1,1,2)-Z(1,1,1));
        Z = cat(3,ones(size(Z,1,2))*z_1_val, Z);
        Z = cat(3,Z,ones(size(Z,1,2))*z_2_val);
        prostate_cube = cat(3,ones(size(prostate_cube,1,2)), prostate_cube);
        prostate_cube = cat(3,prostate_cube, ones(size(prostate_cube,1,2)));
        surface_prostate = isosurface(X, Y, Z, prostate_cube, 0.5);
        p = patch(surface_prostate);
        set(p,'FaceColor',[0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor','none');  
        set(p,'EdgeColor','none');
        camlight;
        lighting gouraud;

        % % add contour
        % info = dicominfo(pathStructureSet);
        % contour = dicomContours(info);
        % try 
        %     contour = deleteContour(contour, 1);
        %     contour = deleteContour(contour, 3); 
        %     contour = deleteContour(contour, 4);
        % catch 
        % end
        % geometricType = "Closed_planar";
        % sliceContour = contour.ROIs(1,:).ContourData{1,1};
        % sliceContour([1:4,6:end]) = [];
        % contour = addContour(contour, 3, "ProstateSlice", sliceContour, geometricType, [0,0,255]);
        % sliceContour = contour.ROIs(2,:).ContourData{1,1};
        % sliceContour([1:4,6:end]) = [];
        % contour = addContour(contour, 4, 'RectumSlice', sliceContour, geometricType, [0,0,255]);
        % contour = deleteContour(contour, 0);
        % contour = deleteContour(contour, 2);
        % plotContour(contour)
        % 
        % % US plane
        % [xi, yi] = meshgrid([-24,24],[-58,-2]);
        % zi = contour.ROIs(1,:).ContourData{1,1}{1}(1,3)*ones(size(xi));
        % ci(:,:,1) = zeros(size(xi));
        % ci(:,:,2) = zeros(size(xi));
        % ci(:,:,3) = ones(size(xi));
        % surf(xi,yi,zi,ci,'EdgeAlpha', 0, 'FaceAlpha', 0.15)

        %plot needles
        n = 1;
        n_idx = cell(1,1);
        s_on_n = [];
        min_z = 0;
        hold on 
        for s = 1:(size(tplan,1)-1)
            if(abs(tplan{s,2}(1)-tplan{s+1,2}(1)) < 1 && abs(tplan{s,2}(2)-tplan{s+1,2}(2)) < 1)
                s_on_n = [s_on_n, s];
            else
                s_on_n = [s_on_n, s];
                n_idx{n} = s_on_n;
                s_on_n = [];
                n = n+1;
            end
            if (tplan{s,2}(3)<min_z)
                min_z = tplan{s,2}(3);
            end
        end
        s_on_n = [s_on_n, s+1];
        n_idx{n} = s_on_n;

        for n = 1:size(n_idx,2)
            ind = n_idx{n};
            point_coordinates = cell2mat(tplan(:,2)');
                needle_x = point_coordinates(1,ind);
            needle_y = point_coordinates(2,ind);
            needle_z = point_coordinates(3,ind);
            if (needle_z(1) < needle_z(end))
                needle_z(1) = min_z - 5;
            elseif (needle_z(1) == needle_z(end))
                needle_x = [needle_x, needle_x];
                needle_y = [needle_y, needle_y];
                needle_z = [needle_z, min_z - 5];
            else
                needle_z(end) = min_z - 5;
            end
            %RGB = squeeze(ind2rgb(max(cindex(ind)),cmap));
            RGB = [0.2,0.2,0.2];
            plot3(needle_x,needle_y,needle_z, 'LineWidth', 1, 'Color', RGB);
        end


        % save figure
        chosen_figure=gcf;
        set(chosen_figure,'PaperUnits','centimeters');
        set(chosen_figure,'PaperSize',[14 14]);
        %saveas(gcf,['..\BRACHYTHERAPY_data\evaluation\3d_visualization_',num2str(id),'_allTargets.pdf'])
        saveas(gcf,['..\BRACHYTHERAPY_data\evaluation\3d_visualization_',num2str(id),'_',targetName,'.pdf'])


        % save gif
        angleRange = -144.5000:2:(360-142.5000);
        gif(['..\BRACHYTHERAPY_data\evaluation\animation_', num2str(id),'_',targetName,'.gif'])
        %gif(['..\BRACHYTHERAPY_data\evaluation\animation_', num2str(id),'_allTargets.gif'])
        for g = 1:length(angleRange)
            % Rotate the camera view
            view(angleRange(g), 30);
            gif
        end

        % planar US view
        img = dicomread(['..\BRACHYTHERAPY_data\',num2str(id-1000),'\IntraOp\IntraOp\MR006.dcm']);
        info = dicominfo(['..\BRACHYTHERAPY_data\',num2str(id-1000),'\IntraOp\IntraOp\MR006.dcm']);
        img = double(img)./255;

        figure
        %imagesc(img)
        warp(X(:,:,slice_num),Y(:,:,slice_num),Z(:,:,slice_num),img);
        colormap(gray)
        hold on 
        plotContour(contour)
        view([0,0,1])
        ax1 = gca;
        axis off

        ax3 = axes;
        hold on
        for s = 1:size(tplan,1)
            point_coordinates = cell2mat(tplan(s,2))';
            if (point_coordinates(3) <= info.ImagePositionPatient(3)+2.5 && point_coordinates(3) >= info.ImagePositionPatient(3)-2.5)
                [Xcirc,Ycirc,Zcirc] = ellipsoid(point_coordinates(1),point_coordinates(2),info.ImagePositionPatient(3),allInverseGradients(s),allInverseGradients(s),0.001);
                x_col = Xcirc-point_coordinates(1);
                y_col = Ycirc-point_coordinates(2);
                r = sqrt(x_col.^2 + y_col.^2)/allInverseGradients(s);
                surf(Xcirc,Ycirc,Zcirc,r,'FaceAlpha',0.4,'EdgeColor','none');
            end
        end      

        warp(X(:,:,slice_num),Y(:,:,slice_num),-50*ones(size(X(:,:,slice_num))),zeros(size(image)));
        view([0,0,1])
        ax3.Color = 'none';
        ax3.View = ax1.View;
        ax3.XLim = ax1.XLim;
        ax3.YLim = ax1.YLim;
        ax3.ZLim = ax1.ZLim;        
        cmap = flipud(hsv);
        cmap(1:156, :) = repmat(cmap(156,:), 156,1);
        colormap(ax3, cmap);
        axis off

         % save figure
        chosen_figure=gcf;
        set(chosen_figure,'PaperUnits','centimeters');
        set(chosen_figure,'PaperSize',[11 8.5]);
        saveas(gcf,['..\BRACHYTHERAPY_data\evaluation\2d_visualization_',num2str(id),'_',num2str(info.ImagePositionPatient(3)),'.pdf'])

end



%% histogram of numbers of sensitive seeds per patient
width = 500;
height = 400;
fontsize = 14;
fontname = 'Times New Roman';
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_numOutliers = [targetName, '_numOutliersPerPatient']; 

    x = evaluation_data{1}.(fieldname_numOutliers);
    edges = [(min(x)-0.5):1:(max(x)+0.5)];
    % Bin the data according to the predefined edges:
    Y = histcounts(x, edges);
    
    binCenters = conv(edges, [0.5, 0.5], 'valid'); % moving average
    [xData, yData] = prepareCurveData( binCenters, Y );
    
    figure('Color', 'w', 'Position', [100, 100, width, height]);
    histogram(x, edges); 
    hold on; 
    grid on;
    xlabel('number of sensitive seeds', 'FontSize', fontsize, 'FontName', fontname)
    ylabel('number of patients', 'FontSize', fontsize, 'FontName', fontname)
    ax = gca; 
    ax.FontSize = fontsize;
    ax.FontName = fontname;
    chosen_figure=gcf;
    set(chosen_figure,'PaperUnits','centimeters');
    set(chosen_figure,'PaperPositionMode','auto');
    set(chosen_figure,'PaperSize',[13 11]);
    set(chosen_figure,'Units','centimeters');
    saveas(gcf,['..\BRACHYTHERAPY_data\evaluation\histogram_',targetName,'.pdf'])
end


%% understand orientation relative to organ surface
% initialize
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_gradPerpendiP = [targetName, '_gradPerpendiP'];
    fieldname_gradPerpendiR = [targetName, '_gradPerpendiR'];
    fieldname_gradPerpendiU = [targetName, '_gradPerpendiU'];
    fieldname_gradParallelP = [targetName, '_gradParallelP'];
    fieldname_gradParallelR = [targetName, '_gradParallelR'];
    fieldname_gradParallelU = [targetName, '_gradParallelU'];
    fieldname_g = [targetName, '_g'];
    evaluation_data{1}.(fieldname_gradPerpendiP) = [];
    evaluation_data{1}.(fieldname_gradPerpendiR) = [];
    evaluation_data{1}.(fieldname_gradPerpendiU) = [];
    evaluation_data{1}.(fieldname_gradParallelP) = [];
    evaluation_data{1}.(fieldname_gradParallelR) = [];
    evaluation_data{1}.(fieldname_gradParallelU) = [];
    evaluation_data{1}.(fieldname_g) = [];
end

% get gradient data
for i = 1:size(patientIds,2)
    for t = 1:length(targetNames)
        targetName = targetNames{t};
        fieldname_gradients = [targetName, '_gradients'];
        fieldname_gradMags = [targetName, '_gradMagnitudes'];
        fieldname_gradPerpendiP = [targetName, '_gradPerpendiP'];
        fieldname_gradPerpendiR = [targetName, '_gradPerpendiR'];
        fieldname_gradPerpendiU = [targetName, '_gradPerpendiU'];
        fieldname_gradParallelP = [targetName, '_gradParallelP'];
        fieldname_gradParallelR = [targetName, '_gradParallelR'];
        fieldname_gradParallelU = [targetName, '_gradParallelU'];
        fieldname_g = [targetName, '_g'];
        fieldname_idx = [targetName, '_outlierIdx'];
        
        gradient_perpendicular_prostate = dot(data{i}.(fieldname_gradients),data{i}.n_prostate, 1).*data{i}.n_prostate;
        gradient_perpendicular_rectum = dot(data{i}.(fieldname_gradients),data{i}.n_rectum, 1).*data{i}.n_rectum;
        gradient_perpendicular_urethra = dot(data{i}.(fieldname_gradients),data{i}.n_urethra, 1).*data{i}.n_urethra;
        gradient_parallel_prostate = data{i}.(fieldname_gradients) - gradient_perpendicular_prostate;
        gradient_parallel_rectum = data{i}.(fieldname_gradients) - gradient_perpendicular_rectum;
        gradient_parallel_urethra = data{i}.(fieldname_gradients) - gradient_perpendicular_urethra;

        % only relevant seeds
        evaluation_data{1}.(fieldname_gradPerpendiP) = [evaluation_data{1}.(fieldname_gradPerpendiP); gradient_perpendicular_prostate(:,data{i}.(fieldname_idx))'];
        evaluation_data{1}.(fieldname_gradPerpendiR) = [evaluation_data{1}.(fieldname_gradPerpendiR); gradient_perpendicular_rectum(:,data{i}.(fieldname_idx))'];
        evaluation_data{1}.(fieldname_gradPerpendiU) = [evaluation_data{1}.(fieldname_gradPerpendiU); gradient_perpendicular_urethra(:,data{i}.(fieldname_idx))'];
        evaluation_data{1}.(fieldname_gradParallelP) = [evaluation_data{1}.(fieldname_gradParallelP); gradient_parallel_prostate(:,data{i}.(fieldname_idx))'];
        evaluation_data{1}.(fieldname_gradParallelR) = [evaluation_data{1}.(fieldname_gradParallelR); gradient_parallel_rectum(:,data{i}.(fieldname_idx))'];
        evaluation_data{1}.(fieldname_gradParallelU) = [evaluation_data{1}.(fieldname_gradParallelU); gradient_parallel_urethra(:,data{i}.(fieldname_idx))'];
        evaluation_data{1}.(fieldname_g) = [evaluation_data{1}.(fieldname_g); repmat({num2str(data{i}.id)}, size(data{i}.(fieldname_idx),2),1)];

    end
end

% plot
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_gradPerpendiP = [targetName, '_gradPerpendiP'];
    fieldname_gradPerpendiR = [targetName, '_gradPerpendiR'];
    fieldname_gradPerpendiU = [targetName, '_gradPerpendiU'];
    fieldname_gradParallelP = [targetName, '_gradParallelP'];
    fieldname_gradParallelR = [targetName, '_gradParallelR'];
    fieldname_gradParallelU = [targetName, '_gradParallelU'];
    fieldname_g = [targetName, '_g'];

    perpendi_norm_p = sqrt(sum(evaluation_data{1, 1}.(fieldname_gradPerpendiP).^2,2));
    perpendi_norm_r = sqrt(sum(evaluation_data{1, 1}.(fieldname_gradPerpendiR).^2,2));
    perpendi_norm_u = sqrt(sum(evaluation_data{1, 1}.(fieldname_gradPerpendiU).^2,2));
    parallel_norm_p = sqrt(sum(evaluation_data{1, 1}.(fieldname_gradParallelP).^2,2));
    parallel_norm_r = sqrt(sum(evaluation_data{1, 1}.(fieldname_gradParallelR).^2,2));
    parallel_norm_u = sqrt(sum(evaluation_data{1, 1}.(fieldname_gradParallelU).^2,2));

    g_ps = [repmat({'P, orthogonal'}, size(evaluation_data{1}.(fieldname_gradPerpendiP),1),1)];
    g_rs = [repmat({'R, orthogonal'}, size(evaluation_data{1}.(fieldname_gradPerpendiR),1),1)];
    g_us = [repmat({'U, orthogonal'}, size(evaluation_data{1}.(fieldname_gradPerpendiU),1),1)];
    g_pp = [repmat({'P, parallel'}, size(evaluation_data{1}.(fieldname_gradParallelP),1),1)];
    g_rp = [repmat({'R, parallel'}, size(evaluation_data{1}.(fieldname_gradParallelR),1),1)];
    g_up = [repmat({'U, parallel'}, size(evaluation_data{1}.(fieldname_gradParallelU),1),1)];

    width = 500;
    height = 400;
    fontsize = 14;
    fontname = 'Times New Roman';
    figure('Color', 'w', 'Position', [100, 100, width, height]);
    boxplot([perpendi_norm_p; parallel_norm_p; perpendi_norm_r;  parallel_norm_r; perpendi_norm_u;  parallel_norm_u], ...
        [g_ps; g_pp; g_rs; g_rp; g_us; g_up]);
    ylabel('gradient magnitude of each component', 'FontSize', fontsize, 'FontName', fontname)
    ax = gca; 
    ax.FontSize = fontsize;
    ax.FontName = fontname;
    chosen_figure=gcf;
    set(chosen_figure,'PaperUnits','centimeters');
    set(chosen_figure,'PaperPositionMode','auto');
    set(chosen_figure,'PaperSize',[13 11]);
    set(chosen_figure,'Units','centimeters');

    % significance
    yt = get(gca, 'YTick');
    axis([xlim 0  max(yt)*1.4])
    xt = get(gca, 'XTick');
    hold on
    [h01,~] = ttest2(perpendi_norm_p, parallel_norm_p,'Alpha',0.01,'Vartype','unequal');
    [h001,~] = ttest2(perpendi_norm_p, parallel_norm_p,'Alpha',0.001,'Vartype','unequal');
    [h0001,~] = ttest2(perpendi_norm_p, parallel_norm_p,'Alpha',0.0001,'Vartype','unequal');
    if (h0001)
        plot(xt([1 1+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([1 1+1])), max(yt)*(1.04+1/10), '*k',  mean(xt([1 1+1]))+0.1, max(yt)*(1.04+1/10), '*k',  mean(xt([1 1+1]))-0.1, max(yt)*(1.04+1/10), '*k')
    elseif (h001)
        plot(xt([1 1+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([1 1+1])), max(yt)*(1.04+1/10), '*k',  mean(xt([1 1+1]))+0.1, max(yt)*(1.04+1/10), '*k')
    elseif (h01)
        plot(xt([1 1+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([1 1+1])), max(yt)*(1.04+1/10), '*k')
    else
        plot(xt([1 1+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([1 1+1])), max(yt)*(1.04+1/10), 'ok')
    end

    [h01,~] = ttest2(perpendi_norm_r, parallel_norm_r,'Alpha',0.01,'Vartype','unequal');
    [h001,~] = ttest2(perpendi_norm_r, parallel_norm_r,'Alpha',0.001,'Vartype','unequal');
    [h0001,~] = ttest2(perpendi_norm_r, parallel_norm_r,'Alpha',0.0001,'Vartype','unequal');
    if (h0001)
        plot(xt([3 3+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([3 3+1])), max(yt)*(1.04+1/10), '*k',  mean(xt([3 3+1]))+0.1, max(yt)*(1.04+1/10), '*k',  mean(xt([3 3+1]))-0.1, max(yt)*(1.04+1/10), '*k')
    elseif (h001)
        plot(xt([3 3+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([3 3+1])), max(yt)*(1.04+1/10), '*k',  mean(xt([3 3+1]))+0.1, max(yt)*(1.04+1/10), '*k')
    elseif (h01)
        plot(xt([3 3+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([3 3+1])), max(yt)*(1.04+1/10), '*k')
    else
        plot(xt([3 3+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([3 3+1])), max(yt)*(1.04+1/10), 'ok')
    end

    [h01,~] = ttest2(perpendi_norm_u, parallel_norm_u,'Alpha',0.01,'Vartype','unequal');
    [h001,~] = ttest2(perpendi_norm_u, parallel_norm_u,'Alpha',0.001,'Vartype','unequal');
    [h0001,~] = ttest2(perpendi_norm_u, parallel_norm_u,'Alpha',0.0001,'Vartype','unequal');
    if (h0001)
        plot(xt([5 5+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([5 5+1])), max(yt)*(1.04+1/10), '*k',  mean(xt([5 5+1]))+0.1, max(yt)*(1.04+1/10), '*k',  mean(xt([5 5+1]))-0.1, max(yt)*(1.04+1/10), '*k')
    elseif (h001)
        plot(xt([5 5+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([5 5+1])), max(yt)*(1.04+1/10), '*k',  mean(xt([5 5+1]))+0.1, max(yt)*(1.04+1/10), '*k')
    elseif (h01)
        plot(xt([5 5+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([5 5+1])), max(yt)*(1.04+1/10), '*k')
    else
        plot(xt([5 5+1]), [1 1]*max(yt)*(1+1/10), '-k',  mean(xt([5 5+1])), max(yt)*(1.04+1/10), 'ok')
    end
    hold off
    saveas(gcf,['..\BRACHYTHERAPY_data\evaluation\gradient_direction_',targetName,'.pdf'])
end


%% Where are the relevant seeds relative to the OARs?
% use distance to mesh to evaluate distance to OAR

% Initialize data
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_dProstate = [targetName, '_dProstate'];
    fieldname_dRectum = [targetName, '_dRectum'];
    fieldname_dUrethra = [targetName, '_dUrethra'];
    evaluation_data{1}.(fieldname_dProstate) = [];
    evaluation_data{1}.(fieldname_dRectum) = [];
    evaluation_data{1}.(fieldname_dUrethra) = [];
end

evaluation_data{1}.dProstateAll = [];
evaluation_data{1}.dRectumAll = [];
evaluation_data{1}.dUrethraAll = [];

for i = 1:size(patientIds,2)

    % get values for all seeds
    evaluation_data{1}.dProstateAll = [evaluation_data{1}.dProstateAll, cell2mat(data{i}.tplan(:,4))'];
    evaluation_data{1}.dRectumAll = [evaluation_data{1}.dRectumAll, cell2mat(data{i}.tplan(:,5))'];
    evaluation_data{1}.dUrethraAll = [evaluation_data{1}.dUrethraAll, cell2mat(data{i}.tplan(:,6))'];
    
    % get value for relevant seeds
    for t = 1:length(targetNames)
        targetName = targetNames{t};
        fieldname_prostate = [targetName, '_dProstate'];
        fieldname_rectum = [targetName, '_dRectum'];
        fieldname_urethra = [targetName, '_dUrethra'];
        fieldname_idx = [targetName, '_outlierIdx'];
        
        distances_prostate =  cell2mat(data{i}.tplan(:,4))'; 
        relevant_distances_prostate = distances_prostate(data{i}.(fieldname_idx));
        evaluation_data{1}.(fieldname_prostate) = [evaluation_data{1}.(fieldname_prostate), relevant_distances_prostate];    

        distances_rectum =  cell2mat(data{i}.tplan(:,5))'; 
        relevant_distances_rectum = distances_rectum(data{i}.(fieldname_idx));
        evaluation_data{1}.(fieldname_rectum) = [evaluation_data{1}.(fieldname_rectum), relevant_distances_rectum];    

        distances_urethra =  cell2mat(data{i}.tplan(:,6))'; 
        relevant_distances_urethra = distances_urethra(data{i}.(fieldname_idx));
        evaluation_data{1}.(fieldname_urethra) = [evaluation_data{1}.(fieldname_urethra), relevant_distances_urethra]; 
    end

end

gP = [repmat({'allSeeds'}, size(evaluation_data{1}.dProstateAll,2),1)];
gR = [repmat({'allSeeds'}, size(evaluation_data{1}.dRectumAll,2),1)];
gU = [repmat({'allSeeds'}, size(evaluation_data{1}.dUrethraAll,2),1)];
dP = [evaluation_data{1}.dProstateAll];
dR = [evaluation_data{1}.dRectumAll];
dU = [evaluation_data{1}.dUrethraAll];

for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_prostate = [targetName, '_dProstate'];
    fieldname_rectum = [targetName, '_dRectum'];
    fieldname_urethra = [targetName, '_dUrethra'];
    gP = [gP; repmat({targetName}, size(evaluation_data{1}.(fieldname_prostate),2),1)];
    gR = [gR; repmat({targetName}, size(evaluation_data{1}.(fieldname_rectum),2),1)];
    gU = [gU; repmat({targetName}, size(evaluation_data{1}.(fieldname_urethra),2),1)];
    dP = [dP, evaluation_data{1}.(fieldname_prostate)];
    dR = [dR, evaluation_data{1}.(fieldname_rectum)];
    dU = [dU, evaluation_data{1}.(fieldname_urethra)];
end

width = 500;
height = 400;
fontsize = 14;
fontname = 'Times New Roman';
figure('Color', 'w', 'Position', [100, 100, width, height]);
boxplot(dP,gP);
ylabel('distance to prostate surface in mm', 'FontSize', fontsize, 'FontName', fontname)
ax = gca; 
ax.FontSize = fontsize;
ax.FontName = fontname;
chosen_figure=gcf;
set(chosen_figure,'PaperUnits','centimeters');
set(chosen_figure,'PaperPositionMode','auto');
set(chosen_figure,'PaperSize',[13 11]);
set(chosen_figure,'Units','centimeters');

yt = get(gca, 'YTick');
axis([xlim 0  max(yt)*1.6])
xt = get(gca, 'XTick');
hold on
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_prostate = [targetName, '_dProstate']; 
    [h01,~] = ttest2(evaluation_data{1}.(fieldname_prostate), evaluation_data{1}.dProstateAll,'Alpha',0.01,'Vartype','unequal');
    [h001,~] = ttest2(evaluation_data{1}.(fieldname_prostate), evaluation_data{1}.dProstateAll,'Alpha',0.001,'Vartype','unequal');
    [h0001,~] = ttest2(evaluation_data{1}.(fieldname_prostate), evaluation_data{1}.dProstateAll,'Alpha',0.0001,'Vartype','unequal');
    if (h0001)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))+0.1, max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))-0.1, max(yt)*(1.04+t/10), '*k')
    elseif (h001)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))+0.1, max(yt)*(1.04+t/10), '*k')
    elseif (h01)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k')
    else
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1  t+1])), max(yt)*(1.04+t/10), 'ok')
    end
end
hold off
saveas(gcf,'..\BRACHYTHERAPY_data\evaluation\seed_position_prostate.pdf')


figure('Color', 'w', 'Position', [100, 100, width, height])
boxplot(dR,gR);
ylabel('distance to rectum surface in mm', 'FontSize', fontsize, 'FontName', fontname)
ax = gca; 
ax.FontSize = fontsize;
ax.FontName = fontname;
chosen_figure=gcf;
set(chosen_figure,'PaperUnits','centimeters');
set(chosen_figure,'PaperPositionMode','auto');
set(chosen_figure,'PaperSize',[13 11]);
set(chosen_figure,'Units','centimeters');
yt = max(evaluation_data{1}.dRectumAll);
axis([xlim 0  max(yt)*1.6])
xt = get(gca, 'XTick');
hold on
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname = [targetName, '_dRectum']; 
    [h01,~] = ttest2(evaluation_data{1}.(fieldname), evaluation_data{1}.dRectumAll,'Alpha',0.01,'Vartype','unequal');
    [h001,~] = ttest2(evaluation_data{1}.(fieldname), evaluation_data{1}.dRectumAll,'Alpha',0.001,'Vartype','unequal');
    [h0001,~] = ttest2(evaluation_data{1}.(fieldname), evaluation_data{1}.dRectumAll,'Alpha',0.0001,'Vartype','unequal');
    if (h0001)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))+0.1, max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))-0.1, max(yt)*(1.04+t/10), '*k')
    elseif (h001)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))+0.1, max(yt)*(1.04+t/10), '*k')
    elseif (h01)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k')
    else
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1  t+1])), max(yt)*(1.04+t/10), 'ok')
    end
end
hold off
saveas(gcf,'..\BRACHYTHERAPY_data\evaluation\seed_position_rectum.pdf')


figure('Color', 'w', 'Position', [100, 100, width, height]);
boxplot(dU,gU);
ylabel('distance to urethra surface in mm', 'FontSize', fontsize, 'FontName', fontname)
ax = gca; 
ax.FontSize = fontsize;
ax.FontName = fontname;
chosen_figure=gcf;
set(chosen_figure,'PaperUnits','centimeters');
set(chosen_figure,'PaperPositionMode','auto');
set(chosen_figure,'PaperSize',[13 11]);
set(chosen_figure,'Units','centimeters');

yt = max(evaluation_data{1}.dUrethraAll);
axis([xlim 0  max(yt)*1.6])
xt = get(gca, 'XTick');
hold on
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname = [targetName, '_dUrethra']; 
    [h01,~] = ttest2(evaluation_data{1}.(fieldname), evaluation_data{1}.dUrethraAll,'Alpha',0.01,'Vartype','unequal');
    [h001,~] = ttest2(evaluation_data{1}.(fieldname), evaluation_data{1}.dUrethraAll,'Alpha',0.001,'Vartype','unequal');
    [h0001,~] = ttest2(evaluation_data{1}.(fieldname), evaluation_data{1}.dUrethraAll,'Alpha',0.0001,'Vartype','unequal');
    if (h0001)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))+0.1, max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))-0.1, max(yt)*(1.04+t/10), '*k')
    elseif (h001)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k',  mean(xt([1 t+1]))+0.1, max(yt)*(1.04+t/10), '*k')
    elseif (h01)
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1 t+1])), max(yt)*(1.04+t/10), '*k')
    else
        plot(xt([1 t+1]), [1 1]*max(yt)*(1+t/10), '-k',  mean(xt([1  t+1])), max(yt)*(1.04+t/10), 'ok')
    end
end
hold off
saveas(gcf,'..\BRACHYTHERAPY_data\evaluation\seed_position_urethra.pdf')


%% Multivariate linear regression
width = 1300;
height = 400;
fontsize = 14;
fontname = 'Times New Roman';
for t = 1:length(targetNames)
    targetName = targetNames{t};
    fieldname_mean = [targetName, '_mean_gradMagnitudes']; % now median
    fieldname_max = [targetName, '_max_gradMagnitudes'];
    fieldname_fval = targetName;

    meanGradMag = [];
    numSeeds = [];
    prostateVol = [];
    variableValue = [];
    activity = [];
    for i = 1:size(patientIds,2)
          meanGradMag = [meanGradMag; data{i}.(fieldname_mean)];
          numSeeds = [numSeeds; size(data{i}.tplan,1)];
          prostateVol = [prostateVol; data{i}.prostateVolume/1000];
          variableValue = [variableValue;data{i}.(fieldname_fval)];
          activity = [activity;data{i}.activity];
    end
    tbl = table(meanGradMag, prostateVol, variableValue, activity);
    mdl = fitlm(tbl,'meanGradMag ~ prostateVol + variableValue + activity','RobustOpts','on')
    figure('Color', 'w', 'Position', [100, 100, width, height]);
    
    % plot(model_mean)
    tiledlayout(1, 3)
    % nexttile
    % plotAdded(mdl, 'numSeeds')
    % title(' ')
    % xlabel('adjusted number of seeds', 'FontSize', fontsize, 'FontName', fontname)
    % ylabel('adjusted median gradient magnitude', 'FontSize', fontsize, 'FontName', fontname)
    % ax = gca; 
    % ax.FontSize = fontsize;
    % ax.FontName = fontname;
    nexttile
    plotAdded(mdl, 'activity')
    title(' ')
    xlabel('adjusted reference air-kerma rate in Gy/s', 'FontSize', fontsize, 'FontName', fontname)
    ylabel('adjusted median gradient magnitude', 'FontSize', fontsize, 'FontName', fontname)
    ax = gca; 
    ax.FontSize = fontsize;
    ax.FontName = fontname;
    nexttile
    plotAdded(mdl, 'prostateVol')
    title(' ')
    xlabel('adjusted prostate volume in cm^3', 'FontSize', fontsize, 'FontName', fontname)
    ylabel('adjusted median gradient magnitude', 'FontSize', fontsize, 'FontName', fontname)
    ax = gca; 
    ax.FontSize = fontsize;
    ax.FontName = fontname;
    nexttile
    plotAdded(mdl, 'variableValue')
    title(' ')

    if (strcmp(targetName,'P_V100'))
        xlabel(['adjusted ',targetName, ' in %'], 'FontSize', fontsize, 'FontName', fontname)
    else
        xlabel(['adjusted ',targetName, ' in Gy'], 'FontSize', fontsize, 'FontName', fontname)
    end
    ylabel('adjusted median gradient magnitude', 'FontSize', fontsize, 'FontName', fontname)
    ax = gca; 
    ax.FontSize = fontsize;
    ax.FontName = fontname;
    chosen_figure=gcf;
    set(chosen_figure,'PaperUnits','centimeters');
    set(chosen_figure,'PaperPositionMode','auto');
    set(chosen_figure,'PaperSize',[32 11]);
    set(chosen_figure,'Units','centimeters');

    saveas(gcf,['..\BRACHYTHERAPY_data\evaluation\regression_',targetName,'.pdf'])

end