% Helper function for mesh_generation. 
% Stacks and pads images for healthy mesh generation.

function imgs = img_merge(folder, format, borderWidth)
    filePattern = fullfile(folder, format);
    bmpFiles = dir(filePattern);
    numFiles = length(bmpFiles);
    addpath(folder);
    
    s = imread(bmpFiles(1).name);
    shape = size(s);
    
    all_data = zeros(shape(1)-(2*borderWidth), shape(2)-(2*borderWidth), numFiles);
    
    for k = 1:numFiles
        Ii = imread(bmpFiles(k).name);
        ext_border = size(Ii);
    
        all_data(:,:,k) = Ii(1+borderWidth:ext_border-borderWidth, 1+borderWidth:ext_border-borderWidth);
    end
    
    imgs = all_data(:,:,1:ext_border-2*(borderWidth));
    
    if max(imgs)>2
        imgs = imgs/255;
    end
end