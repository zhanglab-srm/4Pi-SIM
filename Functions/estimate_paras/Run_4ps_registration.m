function [img_sum,img_B_reg] = Run_4ps_registration(SIM_raw_A, SIM_raw_B)
    % This function performs affine transformation for the 
    % phase-complementary data and return a 6-beam stack.

    disp('  Running registration of rawdata A and B...')
    % the image pair is mirrored
    SIM_raw_B = fliplr(SIM_raw_B);  
    
    % 2D registration
    img_A_sum = zeros(size(SIM_raw_A,1),size(SIM_raw_A,2),1);
    img_B_sum = zeros(size(SIM_raw_A,1),size(SIM_raw_A,2),1);
    for w=1:size(SIM_raw_A,3)
        img_A_sum = img_A_sum + SIM_raw_A(:,:,w);
        img_B_sum = img_B_sum + SIM_raw_B(:,:,w);
    end
    
    % image registration with affine matrix
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 2000;
    fixed = img_A_sum;
    tform = imregtform(img_B_sum,fixed,'affine',optimizer,metric);
    clear img_A_sum img_B_sum
    
    % tranform images
    img_B_reg=zeros(size(SIM_raw_B,1),size(SIM_raw_B,2),size(SIM_raw_B,3));
    for w=1:size(SIM_raw_B,3)
        img_B_reg(:,:,w) = imwarp(SIM_raw_B(:,:,w), tform,'OutputView',imref2d(size(SIM_raw_A)));
    end
    
    % 6-beam stack (simplified spectrum without sidebands, see Supplementary Figure 6)
    img_sum = SIM_raw_A + img_B_reg;

end