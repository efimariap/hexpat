function [SCcomp, SCtheor] = hexpat_1(number_of_meshes, edge_thickness, edge_length, rotation_angles)
    
    % hexpat_1: Function to generate stacked hexagonal meshes moiré pattern.
    % efimariap 2025
    % Inputs:
    %   number_of_meshes - Number of meshes to stack.
    %   edge_thickness   - Thickness of the edges.
    %   edge_length      - Length of the edges.
    %   rotation_angles  - Rotation angles for each mesh.
    %
    % Outputs:
    %   SCcomp   - Computed shade coefficient from the patterns (benchmark values).
    %   SCtheor  - Theoretical shade coefficient.

    % Parameters
    m = number_of_meshes;
    w = edge_thickness;
    d = edge_length;
    theta = rotation_angles;
    
    Lx = 3 * d;
    Ly = ceil(sqrt(3) * d);
  
    
    % Initialize the base hexagonal pattern
    pattern(1:Ly, 1:Lx) = 0;
    pattern(1:w, d/2+1:(d/2+d)) = 1;
    pattern(round(Ly/2-w/2):round(Ly/2+w/2)-1, 2*d-w/2:(2*d+d)) = 1;
    pattern(round(Ly/2), 1:w) = 1;
    pattern(round(Ly/2), (2*d-1):(2*d-1)+w-1) = 1;
    
    % Add diagonal and vertical edges
    for k = 1:d
        % Top-left diagonal
        xt1 = round(k * cos(pi/3));
        yt1 = round(Ly/2 + k * sin(pi/3));
        pattern(yt1, xt1) = 1;
        for iw = 1:100 * w
            xw1 = round(xt1 + 0.01 * iw * sin(pi/6));
            yw1 = round(yt1 - 0.01 * iw * cos(pi/6));
            pattern(yw1, xw1) = 1;
        end
        
        % Bottom-left diagonal
        xt2 = round(k * cos(pi/3));
        yt2 = round(Ly/2 - k * sin(pi/3));
        if yt2 < 1
            yt2 = 1;
        end
        pattern(yt2, xt2) = 1;
        for iw = 1:100 * w
            xw2 = round(xt2 + 0.01 * iw * sin(pi/6));
            yw2 = round(yt2 + 0.01 * iw * cos(pi/6));
            pattern(yw2, xw2) = 1;
        end
        
        % Top-right diagonal
        xt3 = round(Lx/2 + k * (cos(pi/3) - 0.0001));
        yt3 = round(k * sin(pi/3));
        pattern(yt3, xt3) = 1;
        for iw = 1:100 * w
            xw3 = round(xt3 - 0.01 * iw * sin(pi/6));
            yw3 = round(yt3 + 0.01 * iw * cos(pi/6));
            pattern(yw3, xw3) = 1;
        end
        
        % Bottom-right diagonal
        xt4 = round(Lx/2 + k * (cos(pi/3) - 0.000001));
        yt4 = round(Ly - k * sin(pi/3));
        pattern(yt4, xt4) = 1;
        for iw = 1:100 * w
            xw4 = round(xt4 - 0.01 * iw * sin(pi/6));
            yw4 = round(yt4 - 0.01 * iw * cos(pi/6));
            pattern(yw4, xw4) = 1;
        end
    end
    
    % Final adjustments for the pattern
    xt5 = 1;
    yt5 = round(sqrt(3) * (d/2));
    pattern(yt5, xt5) = 1;
    for ir = 1:round(w/2)
        for ic = 1:ceil(ir * tan(pi/6))
            pattern(yt5 + ir, ic) = 1;
            pattern(yt5 - ir, ic) = 1;
        end
    end
    
    % Bottom triangles
    for k = 1:2*w
        pattern(Ly-k+1, round(d/2 + 0.5*k-1):round(d/2 + w - 0.5*k)) = 1;
    end
    for k = 1:2*w
        pattern(Ly-k+1, round((3*d/2) - w + 0.5*k-1):round((3*d/2) - 0.5*k)) = 1;
    end
    
    % Display initial pattern
    impat = mat2gray(pattern);
   % figure, imshow(impat);
    
    % Compute SCcomp and SCtheor for the base pattern
    SCcomp = sum(sum(pattern)) / (Lx * Ly);
    SCtheor = (2/sqrt(3)) * (w/d) - (13/18) * (w/d)^2;
    
    % Create tiled pattern
    nx = 50; ny = 50;
    pattern_total(1:ny*Ly, 1:nx*Lx) = 0;
    for ix = 1:nx
        for iy = 1:ny
            pattern_total(((iy-1)*Ly+1):iy*Ly, ((ix-1)*Lx+1):ix*Lx) = pattern(1:Ly, 1:Lx);
        end
    end
    
    % Process tiled pattern
    impatotal = mat2gray(pattern_total);
    impatotal_mf = medfilt2(impatotal, [1 1]);
    pat_stack(1:ny*Ly, 1:nx*Lx) = impatotal_mf;
    
    % Rotate and stack meshes
    for ip = 1:m-1
        theta1 = theta(ip);
        rotated_pattern = imrotate(impatotal_mf, theta1, 'crop');
        
        % Ensure rotated pattern matches pat_stack dimensions
        if size(rotated_pattern, 1) ~= size(pat_stack, 1) || size(rotated_pattern, 2) ~= size(pat_stack, 2)
            rotated_pattern = imresize(rotated_pattern, [size(pat_stack, 1), size(pat_stack, 2)]);
        end
        
        % Add rotated pattern to stack
        pat_stack(:,:) = pat_stack(:,:) + rotated_pattern;
    end
    
    % Crop the final stacked pattern to 1000x1000
    pat_stack_crop = imcrop(pat_stack, [2000 1500 1500 1500]);
    fig = figure('Visible','off');
    figure, imshow(pat_stack_crop);
   
    
    % Binarize the cropped pattern for SCcomp computation
    pat_stack_crop_bw = imbinarize(pat_stack_crop, 0.1);
    %figure, imshow(pat_stack_crop_bw);
    
    % Compute SCcomp for cropped pattern
    SCcomp = sum(sum(pat_stack_crop_bw)) / (1500 * 1500);
    
    % Compute theoretical SCtheor
    SCtheor = ((2*m)/sqrt(3)) * (w/d) - ((13/18) + (2*(m-1))/(3*sqrt(3))) * m * (w/d)^2;
end


