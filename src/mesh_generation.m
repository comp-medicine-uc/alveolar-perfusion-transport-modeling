%% Current processing pipeline
% Run as Section (in MATLAB) until line 144
% Manually export nodes, elems as .mat files
% Import these files into mesh_processing.py (Python), for each mesh.

%% Array de filenames
clear
folder = '../raw-data/binary_slices/';
format = '*.bmp';


%% Merge images into a NxNxN cubic array
imgs = img_merge(folder, format, 3);


%% Experimental - DICOM image creation
% dicomwritevolume("img.dcm", imgs, 4.5)


%% Image resizing (resolution is 4.5 um per unit in each direction)
tic
RVE_45 = imresize3(imgs, 2, "cubic"); % same and fast
toc
% r2 = imresize3(imgs, 4.5, "linear"); % same
% r3 = imresize3(imgs, 4.5, "box"); % pixelated
% r4 = imresize3(imgs, 4.5, "lanczos2"); % same
% r5 = imresize3(imgs, 4.5, "lanczos3"); % slower


%% Show 
figure
imshow(RVE_45(:,:,1))
title("r")
disp(size(RVE_45))

%% Emphysema-like erosion
RVE_45_adj = zeros(size(RVE_45));
for i=1:size(RVE_45, 3)
    RVE_45_adj(:,:,i) = RVE_45(:,:,i);
end


%% Normal segmentation
offset = 1; %30;
offsetz = 1; %80;
side = 218; %120;

seg_45 = RVE_45_adj(offset:offset+side-1, offset:offset+side-1, offsetz:offsetz+side-1)>25/255;

seg_45_erod1 = zeros(size(seg_45));
seg_45_erod2 = zeros(size(seg_45));
seg_45_erod3 = zeros(size(seg_45));

disk1 = strel("disk", 1);
disk2 = strel("disk", 3);
disk3 = strel("disk", 5);

for i=1:size(seg_45, 3)
    seg_45_erod1(:,:,i) = bwareaopen(imerode(seg_45(:,:,i), disk1), 4, 4);
    seg_45_erod2(:,:,i) = bwareaopen(imerode(imerode(seg_45(:,:,i), disk1), disk1), 8, 4);
    seg_45_erod3(:,:,i) = bwareaopen(imerode(imerode(imerode(seg_45(:,:,i), disk1), disk1), disk1), 8, 4);
end

seg_45_erod1 = logical(seg_45_erod1);
seg_45_erod2 = logical(seg_45_erod2);
seg_45_erod3 = logical(seg_45_erod3);


%% Porosity calculation and plotting
porosity_vol_normal = 1 - nnz(seg_45)/numel(seg_45);
porosity_vol_1 = 1 - nnz(seg_45_erod1)/numel(seg_45_erod1);
porosity_vol_2 = 1 - nnz(seg_45_erod2)/numel(seg_45_erod2);
porosity_vol_3 = 1 - nnz(seg_45_erod3)/numel(seg_45_erod3);

disp("Normal porosity = " + num2str(porosity_vol_normal))
disp("Emphysema I porosity = " + num2str(porosity_vol_1))
disp("Emphysema II porosity = " + num2str(porosity_vol_2))
disp("Emphysema III porosity = " + num2str(porosity_vol_3))

layer = 40;

figure
montage([RVE_45(:,:,layer), RVE_45_adj(:,:,layer)])
figure
montage([seg_45(:,:,layer), seg_45_erod1(:,:,layer), seg_45_erod2(:,:,layer), seg_45_erod3(:,:,layer)])


%% Plotting
}
for i=1:3
    figure
    % montage([RVE_45(:,:,2*(i-1)+1), seg_45(:,:,2*(i-1)+1)])
    imshow(seg_45_erod3(:,:,2*(i-1)+1))
end


%% 3D meshing with iso2mesh
opt.radbound = 6;
[nodes, elems] = v2s(seg_45, 1, opt, 'cgalmesh');
[nodes, elems] = meshcheckrepair(nodes, elems, 'meshfix');
nodes = sms(nodes, elems, 10, 0.8, 'lowpass');
[nodes, elems] = meshcheckrepair(nodes, elems, 'meshfix');

opt.radbound = 6;
[nodes1, elems1] = v2s(seg_45_erod1, 1, opt, 'cgalmesh');
[nodes1, elems1] = meshcheckrepair(nodes1, elems1, 'meshfix');
nodes1 = sms(nodes1, elems1, 10, 0.8, 'lowpass');
[nodes1, elems1] = meshcheckrepair(nodes1, elems1, 'meshfix');

opt.radbound = 6;
[nodes2, elems2] = v2s(seg_45_erod2, 1, opt, 'cgalmesh');
[nodes2, elems2] = meshcheckrepair(nodes2, elems2, 'meshfix');
nodes2 = sms(nodes2, elems2, 10, 0.8, 'lowpass');
[nodes2, elems2] = meshcheckrepair(nodes2, elems2, 'meshfix');

opt.radbound = 6;
[nodes3, elems3] = v2s(seg_45_erod3, 1, opt, 'cgalmesh');
[nodes3, elems3] = meshcheckrepair(nodes3, elems3, 'meshfix');
nodes3 = sms(nodes3, elems3, 10, 0.8, 'lowpass');
[nodes3, elems3] = meshcheckrepair(nodes3, elems3, 'meshfix');


%% Expansion of nodes to match image scale 
% (each voxel is a square of 2.25 um of side length)

nodes = nodes * 2.25;
nodes1 = nodes1 * 2.25;
nodes2 = nodes2 * 2.25;
nodes3 = nodes3 * 2.25;


%% Emphysematous RVE 3D plotting

figure
title("Normal lung RVE mesh")
plotmesh(nodes, elems)
figure
title("Emphysema I lung RVE mesh")
plotmesh(nodes1, elems1)
figure
title("Emphysema II lung RVE mesh")
plotmesh(nodes2, elems2)
figure
title("Emphysema III lung RVE mesh")
plotmesh(nodes3, elems3)


%% Slice-by-slice erosion and dilation (experimental)

% side = 2.25*100;
% offset = 29;
% offsetz = 79;
% layers = 1:side;
% alpha = 0.2;
% threshs = [25];

% r_eroded_1 = zeros([side,side,side]);
% r_eroded_2 = zeros([side,side,side]);

% seg = r(offset:side+offset-1,offset:side+offset-1,offsetz:side+offsetz-1)>26/255;

% radius1 = 5;
% radius2 = 10;
% disk1 = strel("disk", radius1);
% disk2 = strel("disk", radius2);

% for i=1:size(seg, 3)

%     temp = seg(:,:,i);
%     r_skeleton = imdilate(bwmorph(temp,"skel",Inf), strel("disk", 1));
%     r_eroded_1(:,:,i) = bwmorph(bitor(imerode(temp, disk1), r_skeleton), "clean");
%     r_eroded_2(:,:,i) = bwmorph(bitor(imerode(temp, disk2), r_skeleton), "clean");
    
%     % Sanity check
%     % 
%     % if i==30
%     %     figure
%     %     out = seg(:,:,i)*alpha + r(offset:side+offset,offset:side+offset,offset+35+i-1)*1.5;
%     %     imshow(out)
%     %     title("normal")
%     %     figure
%     %     out = r_eroded_1(:,:,i)*alpha + r(offset:side+offset,offset:side+offset,offset+35+i-1)*1.5;
%     %     imshow(out)
%     %     title("radius" + radius1)
%     %     figure
%     %     out = r_eroded_2(:,:,i)*alpha + r(offset:side+offset,offset:side+offset,offset+35+i-1)*1.5;
%     %     imshow(out)
%     %     title("radius" + radius2)
%     %     figure
%     %     out = r_skeleton*alpha + r(offset:side+offset,offset:side+offset,offset+35+i-1)*1.5;
%     %     imshow(out)
%     % end

% end

% seg = imresize3(seg, 1/2.25, "nearest");
% r_eroded_1 = imresize3(r_eroded_1, 1/2.25, "nearest");
% r_eroded_2 = imresize3(r_eroded_2, 1/2.25, "nearest");
% r_skeleton = imresize3(r_skeleton, 1/2.25, "nearest"); % Skeletonized image (experimental)

% disp("Normal volume porosity = " + num2str(1 - nnz(seg)/numel(seg)));
% disp("Dilated 1 volume porosity = " + num2str(1 - nnz(r_eroded_1)/numel(r_eroded_1)))
% disp("Dilated 2 volume porosity = " + num2str(1 - nnz(r_eroded_2)/numel(r_eroded_2)))

% j = 30;

% bg = (1/1.5)*imresize3(r(offset:side+offset,offset:side+offset,offsetz:side+offsetz-1), 1/2.26, "cubic");

% figure
% out = seg(:,:,j)*alpha + bg(:,:,j)*1.5;
% imshow(out)
% title("normal")
% figure
% out = r_eroded_1(:,:,j)*alpha + bg(:,:,j)*1.5;
% imshow(out)
% title("radius" + radius1)
% figure
% out = r_eroded_2(:,:,j)*alpha + bg(:,:,j)*1.5;
% imshow(out)
% title("radius" + radius2)
% figure
% out = r_skeleton*alpha + bg(:,:,100)*1.5;
% imshow(out)

%% Outputs of previous (commented) block
% Normal volume porosity = 0.64991
% Dilated 1 volume porosity = 0.77608
% Dilated 2 volume porosity = 0.85091


%% Data binarization loop (experimental)
side = 100;
offset = 13;
offsetz = 35;
layers = 1:side;
alpha = 0.2;
threshs = [25];
ps = zeros(3, length(layers));
pvols = zeros(length(threshs),1);

r_normal =  imresize3(r, 1/2.25, "cubic");
r_eroded_a =  imerode(r, strel('sphere',3));
r_eroded_b =  imerode(r, strel('sphere',4));
% r_dilated = imdilate(r, strel('sphere',3));

x0=0.3*1920;
y0=0.6*1080/2;
width=700;
height=700;

% for i=1:length(threshs)
%     seg_dilated = r_dilated > threshs(i)/255;
%     seg_dilated = seg_dilated(offset:side+offset,offset:side+offset,offset+35:side+offset+35);
%     for j=1:length(layers)
%         porosityarea = 1-nnz(seg_dilated(:,:,layers(j)))/numel(seg_dilated(:,:,layers(j)));
%         ps(1,j) = porosityarea;
%         if j==1
%             porosityvol = 1-nnz(seg_dilated)/numel(seg_dilated);
%             figure
%             out = seg_dilated(:,:,layers(j))*alpha + r(offset:side+offset,offset:side+offset,layers(j+offset+35))*1.5;
%             imshow(out)
%             title("threshold = " + num2str(threshs(i)) + ", porosityvol = "+ num2str(porosityvol)+ ", porosityarea = "+ num2str(porosityarea))
%             set(gcf,'position',[x0,y0,width,height])
%         end
%     end
%     pvols(1) = porosityvol;
% end

for i=1:length(threshs)
    seg_normal = r_normal > threshs(i)/255;
    seg_normal = seg_normal(offset:side+offset,offset:side+offset,offset+35:side+offset+35);
    for j=1:length(layers)
        porosityarea = 1-nnz(seg_normal(:,:,layers(j)))/numel(seg_normal(:,:,layers(j)));
        ps(1,j) = porosityarea;
        if j==1
            porosityvol = 1-nnz(seg_normal)/numel(seg_normal);
            figure
            out = seg_normal(:,:,layers(j))*alpha + r(offset:side+offset,offset:side+offset,layers(j+offset+35))*1.5;
            imshow(out)
            title("threshold = " + num2str(threshs(i)) + ", porosityvol = "+ num2str(porosityvol)+ ", porosityarea = "+ num2str(porosityarea))
            set(gcf,'position',[x0,y0,width,height])
        end
    end
    pvols(1) = porosityvol;
end

for i=1:length(threshs)
    seg_eroded_a = r_eroded_a > threshs(i)/255;
    seg_eroded_a = seg_eroded_a(offset:side+offset,offset:side+offset,offset+35:side+offset+35);
    for j=1:length(layers)
        % seg_eroded(:,:,j) = bwareaopen(seg_eroded(:,:,j), 8, 4);
        porosityarea = 1-nnz(seg_eroded_a(:,:,layers(j)))/numel(seg_eroded_a(:,:,layers(j)));
        ps(2,j) = porosityarea;
        if j==1
            porosityvol = 1-nnz(seg_eroded_a)/numel(seg_eroded_a);
            figure
            out = seg_eroded_a(:,:,layers(j))*alpha + r(offset:side+offset,offset:side+offset,layers(j+offset+35))*1.5;
            imshow(out)
            title("threshold = " + num2str(threshs(i)) + ", porosityvol = "+ num2str(porosityvol)+ ", porosityarea = "+ num2str(porosityarea))
            set(gcf,'position',[x0,y0,width,height])
        end
    end
    pvols(2) = porosityvol;
end

for i=1:length(threshs)
    seg_eroded_b = r_eroded_b > threshs(i)/255;
    seg_eroded_b = seg_eroded_b(offset:side+offset,offset:side+offset,offset+35:side+offset+35);
    for j=1:length(layers)
        % seg_eroded(:,:,j) = bwareaopen(seg_eroded(:,:,j), 8, 4);
        porosityarea = 1-nnz(seg_eroded_b(:,:,layers(j)))/numel(seg_eroded_b(:,:,layers(j)));
        ps(3,j) = porosityarea;
        if j==1
            porosityvol = 1-nnz(seg_eroded_b)/numel(seg_eroded_b);
            figure
            out = seg_eroded_b(:,:,layers(j))*alpha + r(offset:side+offset,offset:side+offset,layers(j+offset+35))*1.5;
            imshow(out)
            title("threshold = " + num2str(threshs(i)) + ", porosityvol = "+ num2str(porosityvol)+ ", porosityarea = "+ num2str(porosityarea))
            set(gcf,'position',[x0,y0,width,height])
        end
    end
    pvols(3) = porosityvol;
end

% threshold = 25 gives porosity = 0.65 globally, which agrees with Sarabia &
% Hurtado. threshold = 25 was chosen for a sub-cube, even though it gives
% 0.668 porosity.

figure
hold on
for i=1:3
    plot(1:side, ps(i,:), 'DisplayName', "threshold = " + num2str(threshs(1)) + ", avg porosity = " + num2str(pvols(i)))
end
legend
hold off


%% RVE tesellation (deprecated)

% % Lower layer is indexed as follows (as seen from positive z axis)
% % 1 2
% % 3 4
% lower = zeros(2*side,2*side,side);
% for i = 1:side
%     layer1 = seg(1:side,1:side,i);
%     layer2 = fliplr(layer1);
%     layer3 = flipud(layer1);
%     layer4 = rot90(layer1,2);
%     lower(1:side,1:side,i) = layer1;
%     lower(side+1:2*side,1:side,i) = layer3(:,:);
%     lower(1:side,side+1:2*side,i) = layer2(:,:);
%     lower(side+1:2*side,side+1:2*side,i) = layer4(:,:);
% end
% 
% full = zeros(2*side,2*side,2*side);
% for j = 1:side
%     full(:,:,j) = lower(:,:,j);
% end
% for k = side+1:2*side
%     full(:,:,k) = lower(:,:,2*side-k+1);
% end


%% Single-RVE mesh generation and repairing (works!)
% opt.radbound = 6;
% [nodes, elems] = v2s(seg, 1, opt, 'cgalmesh');
% [nodes, elems] = meshcheckrepair(nodes, elems, 'meshfix');
% nodes = sms(nodes, elems, 10, 0.8, 'lowpass');
% [nodes, elems] = meshcheckrepair(nodes, elems, 'meshfix');


%% Multi-meshing (experimental)
% disp('Dilated mesh generation')
% opt.radbound = 6;
% [nodes_d, elems_d] = v2s(seg_dilated, 1, opt, 'cgalmesh');
% [nodes_d, elems_d] = meshcheckrepair(nodes_d, elems_d, 'meshfix');
% nodes_d = sms(nodes_d, elems_d, 10, 0.8, 'lowpass');
% [nodes_d, elems_d] = meshcheckrepair(nodes_d, elems_d, 'meshfix');
% %
% disp('Normal mesh generation')
% opt.radbound = 6;
% [nodes_n, elems_n] = v2s(seg_normal, 1, opt, 'cgalmesh');
% [nodes_n, elems_n] = meshcheckrepair(nodes_n, elems_n, 'meshfix');
% nodes_n = sms(nodes_n, elems_n, 10, 0.8, 'lowpass');
% [nodes_n, elems_n] = meshcheckrepair(nodes_n, elems_n, 'meshfix');

% disp('Eroded mesh generation')
% opt.radbound = 6;
% [nodes_e, elems_e] = v2s(r_eroded, 1, opt, 'cgalmesh');
% [nodes_e, elems_e] = meshcheckrepair(nodes_e, elems_e, 'meshfix');
% nodes_e = sms(nodes_e, elems_e, 10, 0.8, 'lowpass');
% [nodes_e, elems_e] = meshcheckrepair(nodes_e, elems_e, 'meshfix');

disp('Eroded 1 mesh generation')
opt.radbound = 6;
[nodes_e_1, elems_e_1] = v2s(r_eroded_1, 1, opt, 'cgalmesh');
[nodes_e_1, elems_e_1] = meshcheckrepair(nodes_e_1, elems_e_1, 'meshfix');
nodes_e_1 = sms(nodes_e_1, elems_e_1, 10, 0.8, 'lowpass');
[nodes_e_1, elems_e_1] = meshcheckrepair(nodes_e_1, elems_e_1, 'meshfix');

disp('Eroded 2 mesh generation')
opt.radbound = 6;
[nodes_e_2, elems_e_2] = v2s(r_eroded_2, 1, opt, 'cgalmesh');
[nodes_e_2, elems_e_2] = meshcheckrepair(nodes_e_2, elems_e_2, 'meshfix');
nodes_e_2 = sms(nodes_e_2, elems_e_2, 10, 0.8, 'lowpass');
[nodes_e_2, elems_e_2] = meshcheckrepair(nodes_e_2, elems_e_2, 'meshfix');


%% Mesh final node resizing
nodes = nodes * 2.25;
% nodes_d = nodes_d * 2.25;
% nodes_n = nodes_n * 2.25;
% nodes_e = nodes_e * 2.25;

nodes_e_1 = nodes_e_1 * 2.25;
nodes_e_2 = nodes_e_2 * 2.25;

%% Emphysematous RVE 3D plotting
figure
title("Normal lung RVE mesh")
plotmesh(nodes, elems)
figure
title("Emphysema I lung RVE mesh")
plotmesh(nodes_e_1, elems_e_1)
figure
title("Emphysema II lung RVE mesh")
plotmesh(nodes_e_2, elems_e_2)
% disp(["nodes length = ", len(nodes)]);


%% MATLAB mesh face flattening
% Inlet and outlet face flattening reduces numerical errors
% in both boundaries. 

mod_elems = elems3;
mod_nodes = nodes3;

min_x = min(mod_nodes(:,1));
disp(["min_x = ", min_x]);
min_y = min(mod_nodes(:,2));
disp(["min_y = ", min_y]);
min_z = min(mod_nodes(:,3));
disp(["min_z = ", min_z]);
max_x = max(mod_nodes(:,1));
disp(["max_x = ", max_x]);
max_y = max(mod_nodes(:,2));
disp(["max_y = ", max_y]);
max_z = max(mod_nodes(:,3));
disp(["max_z = ", max_z]);


%% Flattening step
tol = 2;
d = 2.5;
for i = 1:length(mod_nodes)
    if abs(mod_nodes(i,1)-min_x)<tol && abs(mod_nodes(i,2)-min_y)>2*tol && abs(mod_nodes(i,2)-max_y)>2*tol && abs(mod_nodes(i,3)-min_z)>2*tol && abs(mod_nodes(i,3)-max_z)>2*tol
        mod_nodes(i,1)=min_x-d;
    elseif abs(mod_nodes(i,1)-max_x)<tol && abs(mod_nodes(i,2)-min_y)>2*tol && abs(mod_nodes(i,2)-max_y)>2*tol && abs(mod_nodes(i,3)-min_z)>2*tol && abs(mod_nodes(i,3)-max_z)>2*tol
        mod_nodes(i,1)=max_x+d;
    % elseif mod_nodes(i,2)<tol && ~(mod_nodes(i,1)<2*tol||mod_nodes(i,1)>len-2*tol) && ~(mod_nodes(i,3)<2*tol||mod_nodes(i,3)>len-2*tol)
    %     mod_nodes(i,2)=min_y-d;
    % elseif mod_nodes(i,2)>len-tol && ~(mod_nodes(i,1)<2*tol||mod_nodes(i,1)>len-2*tol) && ~(mod_nodes(i,3)<2*tol||mod_nodes(i,3)>len-2*tol)
    %     mod_nodes(i,2)=max_y+d;
    % elseif mod_nodes(i,3)<tol && ~(mod_nodes(i,1)<tol||mod_nodes(i,1)>len-tol) && ~(mod_nodes(i,2)<tol||mod_nodes(i,2)>len-tol)
    %     mod_nodes(i,3)=min_z-d;
    % elseif mod_nodes(i,3)>len-tol && ~(mod_nodes(i,1)<tol||mod_nodes(i,1)>len-tol) && ~(mod_nodes(i,2)<tol||mod_nodes(i,2)>len-tol)
    %     mod_nodes(i,3)=max_z+d;
    end
end


%% Plot after flattening
figure
title("Mesh after flattening")
plotmesh(mod_nodes, mod_elems)