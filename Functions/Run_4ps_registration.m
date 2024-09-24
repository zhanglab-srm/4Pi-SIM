function [img_sum,img_right_reg] = Run_4ps_registration(SIM_raw_A, SIM_raw_B)

SIM_raw_B=fliplr(SIM_raw_B);

%sum of left and right image
img_left_sum=zeros(size(SIM_raw_A,1),size(SIM_raw_A,2),1);
img_right_sum=zeros(size(SIM_raw_A,1),size(SIM_raw_A,2),1);
for w=1:size(SIM_raw_A,3)
    img_left_sum=img_left_sum+SIM_raw_A(:,:,w);
    img_right_sum=img_right_sum+SIM_raw_B(:,:,w);
end

% image registration and make tform
[optimizer,metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 2000;
fixed = img_left_sum;
tform = imregtform(img_right_sum,fixed,'affine',optimizer,metric);
clear img_left_sum img_right_sum

% Tranform images
img_right_reg=zeros(size(SIM_raw_B,1),size(SIM_raw_B,2),size(SIM_raw_B,3));
for w=1:size(SIM_raw_B,3)
    img_right_reg(:,:,w) = imwarp(SIM_raw_B(:,:,w), tform,'OutputView',imref2d(size(SIM_raw_A)));
end
img_sum=SIM_raw_A + img_right_reg;

end