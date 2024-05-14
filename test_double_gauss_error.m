% clear;
% clc;
% CoreNum = 6; %设定机器CPU核心数量
% if isempty(gcp('nocreate')) %如果并行未开启
%     p=parpool(CoreNum);
% end

tic;
%折射率：587.6nm单波长下
sk2=1.6073788586;
sk16=1.6204079331;
f5=1.6034171782;
%折射率：555.0nm单波长下
bk7 = 1.5185223876;

air=1;
opt = opt_model();
a = opt.seq_model;
a.source_thick = 300;% inf=1e10
a.source_field_objheight = [0,0];
a.setEntrancePupil(33.33); %设置入瞳直径

%% 双高斯透镜
a.initialInterface()
a.addInterface('Conic', 54.153, 8.747, sk2, 1);
a.addInterface('Conic', 152.522, 0.5, air, 1);
a.addInterface('Conic', 35.951, 14, sk16, 1);
a.addInterface('Conic', inf, 3.777, f5, 1);
a.addInterface('Conic', 22.270, 14.253, air, 1);
a.addInterface('Conic', inf, 12.428, air, 1);
a.set_stop();
a.addInterface('Conic', -25.685, 3.777, f5, 1);
a.addInterface('Conic', inf, 10.834, sk16, 1);
a.addInterface('Conic', -36.980, 0.5, air, 1);
a.addInterface('Conic', 196.417, 6.858, sk16, 1);
a.addInterface('Conic', -67.148, 94.071, air, 1);
%% 像面
img_source = imread('Demo_picture_640x480.bmp');
img_source_height = 150;

img_plane_sample = 32;
psf_x_num = 4;
psf_y_num = 5;

% 计算当前设定的视场高度下，光源图片对应的x,y物理尺寸a.source_xwidth a.source_ywidth
img_row_pix = size(img_source,1);
img_col_pix = size(img_source,2);
img_color_num = size(img_source,3);
a.source_pixel_xnum = img_col_pix;
a.source_pixel_ynum = img_row_pix;
a.source_xwidth = img_source_height*img_col_pix/img_row_pix;
a.source_ywidth = img_source_height;

% 计算设定的psf点数下，不同psf点源之间的x,y物理距离delta_fieldheight_x delta_fieldheight_y
psf_tol_num = psf_x_num*psf_y_num;
chiefray_obj = zeros(psf_tol_num,3);
chiefray_img = zeros(psf_tol_num,3);
delta_fieldheight_x = a.source_xwidth/psf_x_num;
delta_fieldheight_y = a.source_ywidth/psf_y_num;

% 计算像面位图像元的默认物理尺寸img_pix_size，单个psf图像的默认物理尺寸img_x_width img_y_width
a.field_height_x = a.source_field_objheight(1);
a.field_height_y = a.source_field_objheight(2);
a.update();
img_field_objheight = a.chiefray_s_back;
% 默认情况下，像面只需给出采样网格数，网格大小的计算如下
x0 = a.p_s_back(1,1);y0 = a.p_s_back(1,2);
x1 = a.p_s_back(2,1);y1 = a.p_s_back(2,2);
img_pix_size = sqrt(2*((x1-x0)^2 + (y1-y0)^2))/2;%默认像元尺寸
img_sample_xnum = img_plane_sample;
img_sample_ynum = img_plane_sample;
img_x_width = img_sample_xnum*img_pix_size;%单个PSF图像的尺寸
img_x_halfwidth = img_x_width/2;
img_y_width = img_sample_ynum*img_pix_size;
img_y_halfwidth = img_y_width/2;
img_tol_xwidth = img_pix_size*img_col_pix;%整个像面图像的尺寸
img_tol_ywidth = img_pix_size*img_row_pix;

% 计算所有psf图像
ipsf = 0;
I_ipsf = zeros(img_sample_ynum,img_sample_xnum,psf_tol_num);
for ifieldy = 1:psf_y_num
    for ifieldx = 1:psf_x_num
        ipsf = ipsf + 1;
        a.field_height_x = a.source_field_objheight(1) - a.source_xwidth/2 - delta_fieldheight_x/2 + ifieldx*delta_fieldheight_x;
        a.field_height_y = a.source_field_objheight(2) + a.source_ywidth/2 + delta_fieldheight_y/2 - ifieldy*delta_fieldheight_y;
        a.update();
        chiefray_obj(ipsf,:) = a.chiefray_s_start;
        chiefray_img(ipsf,:) = a.chiefray_s_back;
        x_min = chiefray_img(ipsf,1)-img_x_halfwidth;x_max = chiefray_img(ipsf,1)+img_x_halfwidth;
        y_min = chiefray_img(ipsf,2)-img_y_halfwidth;y_max = chiefray_img(ipsf,2)+img_y_halfwidth;
        I_ipsf(:,:,ipsf) = I_add(a,x_min,x_max,img_sample_xnum,y_min,y_max,img_sample_ynum);
    end
end

% 计算光源位图中，每个采样视场点源所在的像元索引
x_detect_spot = linspace(-a.source_xwidth/2,a.source_xwidth/2,a.source_pixel_xnum);
y_detect_spot = linspace(a.source_ywidth/2,-a.source_ywidth/2,a.source_pixel_ynum);%注意第一行对应+ymax
index_cheifray_obj = zeros(psf_tol_num,2);
for icheifray_obj = 1:psf_tol_num
    [~,x_near] = min(abs(x_detect_spot(:) - chiefray_obj(icheifray_obj,1) + a.source_field_objheight(1)));
    [~,y_near] = min(abs(y_detect_spot(:) - chiefray_obj(icheifray_obj,2) + a.source_field_objheight(2)));
    index_cheifray_obj(icheifray_obj,:)  = [y_near,x_near];%注意：x对应列，y对应行
end

% 绘制psf网格图，网格图中显示的单个psf图的中心位置和点源的位置一致
I_psf_grid = zeros(img_row_pix,img_col_pix);
for ipsf = 1:psf_tol_num
    I_psf_grid_temp = zeros(img_row_pix,img_col_pix);
    x_left = index_cheifray_obj(ipsf,2) - floor(img_plane_sample/2);
    x_right = x_left + img_plane_sample - 1;
    y_left = index_cheifray_obj(ipsf,1) - floor(img_plane_sample/2);
    y_right = y_left + img_plane_sample - 1;
    I_psf_grid_temp(y_left:y_right,x_left:x_right) = I_ipsf(:,:,ipsf);
    I_psf_grid = I_psf_grid + I_psf_grid_temp;
end
figure(1);
imagesc(I_psf_grid);
title('PSF网格');
axis equal tight;
colorbar;

img_afterSINC = zeros(img_row_pix,img_col_pix,img_color_num);
img_afterSINC_unflip = zeros(img_row_pix,img_col_pix,img_color_num);
for icolor = 1:img_color_num
    img_source_color = double(img_source(:,:,icolor));
    %% 先把整张图和不同视场的psf核进行卷积
    img_convol = zeros(img_row_pix,img_col_pix,psf_tol_num);
    for ifield = 1:psf_tol_num
        img_convol(:,:,ifield) = conv2(img_source_color,I_ipsf(:,:,ifield),'same');
    end

    %% 直接计算整张图每个像素的插值
    img_afterSINC_color_unflip = zeros(img_row_pix,img_col_pix);
    img_xhw = img_tol_xwidth/2; 
    img_yhw = img_tol_ywidth/2;
    % 计算像面位图中，每个采样视场点源对应像面所在的像元索引
    x_detect_spot = linspace(-img_xhw,img_xhw,img_col_pix);
    y_detect_spot = linspace(img_yhw,-img_yhw,img_row_pix);%注意第一行对应+ymax
    index_cheifray_img = zeros(psf_tol_num,2);
    for icheifray_img = 1:psf_tol_num
        [~,x_near] = min(abs(x_detect_spot(:) - chiefray_img(icheifray_img,1) + img_field_objheight(1)));
        [~,y_near] = min(abs(y_detect_spot(:) - chiefray_img(icheifray_img,2) + img_field_objheight(2)));
        index_cheifray_img(icheifray_img,:)  = [y_near,x_near];%注意：x对应列，y对应行
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 并行计算
    y = 1:img_row_pix;
    x = 1:img_col_pix;
    [X,Y] = meshgrid(x,y);
    source_yhw_patch = a.source_ywidth/2/psf_y_num;
    source_xhw_patch = a.source_xwidth/2/psf_x_num;
    source_pix_size = a.source_ywidth/img_row_pix;
    for ii = 1:img_row_pix*img_col_pix
            S_sum = 0;
            img_afterSINC_color_p = 0;
            i_pix = Y(ii);
            j_pix = X(ii);
            for i_psf = 1:psf_tol_num
                y_s = source_pix_size*(i_pix - index_cheifray_obj(i_psf,1))/source_yhw_patch/2; %注意这里要多除以一个2
                x_s = source_pix_size*(j_pix - index_cheifray_obj(i_psf,2))/source_xhw_patch/2;
                S_pix = sinc(x_s)*sinc(y_s);
                img_afterSINC_color_p = img_afterSINC_color_p + S_pix*img_convol(i_pix,j_pix,i_psf);
                S_sum = S_sum + S_pix;
            end
            img_afterSINC_color_unflip(ii) = img_afterSINC_color_p/S_sum;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if a.opt_model.ParaxialModel.paraxial_magnification<0
        img_afterSINC_color1 = flipud(img_afterSINC_color_unflip);
        img_afterSINC_color = fliplr(img_afterSINC_color1);
    end
    img_afterSINC_unflip(:,:,icolor) = img_afterSINC_color_unflip;
    img_afterSINC(:,:,icolor) = img_afterSINC_color;
end

img_afterSINC = uint8(img_afterSINC);
img_afterSINC_unflip_uint = uint8(img_afterSINC_unflip);

if a.opt_model.ParaxialModel.paraxial_magnification<0
    img_source1 = flipud(img_source);
    img_source_for_vs = fliplr(img_source1);
end

%% 畸变处理
xmin_distort_obj = a.source_field_objheight(1) - a.source_xwidth/2;
xmax_distort_obj = a.source_field_objheight(1) + a.source_xwidth/2;
ymin_distort_obj = a.source_field_objheight(2) - a.source_ywidth/2;
ymax_distort_obj = a.source_field_objheight(2) + a.source_ywidth/2;
% 计算畸变网格点
x_detect_spot = linspace(xmin_distort_obj,xmax_distort_obj,a.source_pixel_xnum);
y_detect_spot = linspace(ymax_distort_obj,ymin_distort_obj,a.source_pixel_ynum);%注意第一行对应+ymax
[Y_detect_spot,X_detect_spot] = meshgrid(y_detect_spot,x_detect_spot);
distortray_obj = zeros(img_row_pix*img_col_pix,3);
for idistortray = 1:img_row_pix*img_col_pix
    distortray_obj(idistortray,:) = [X_detect_spot(idistortray),Y_detect_spot(idistortray),0];
end
a.distortarray_s_start = distortray_obj;
a.update();
distortray_img = a.distortarray_s_back;
% 畸变网格点与模糊图像灰度值对应
c1 = reshape((img_afterSINC_unflip(:,:,1))',[img_row_pix*img_col_pix,1]);
c2 = reshape((img_afterSINC_unflip(:,:,2))',[img_row_pix*img_col_pix,1]);
c3 = reshape((img_afterSINC_unflip(:,:,3))',[img_row_pix*img_col_pix,1]);
distortray_img_3c = [distortray_img(:,1),distortray_img(:,2),c1,c2,c3];
% 对畸变网格点进行像元归属计算
x_min = img_field_objheight(1) - img_tol_xwidth/2;
x_max = img_field_objheight(1) + img_tol_xwidth/2;
y_min = img_field_objheight(2) - img_tol_ywidth/2;
y_max = img_field_objheight(2) + img_tol_ywidth/2;
x_pixels = img_col_pix;
y_pixels = img_row_pix;
s_img = distortray_img_3c;
nray = size(s_img,1);
x_detect_spot = linspace(x_min,x_max,x_pixels);
y_detect_spot = linspace(y_max,y_min,y_pixels);
img_distort_temp = zeros(img_row_pix,img_col_pix,3);
pixel_ray_count = zeros(img_row_pix,img_col_pix,3);
for j = 1:nray
    if(s_img(j,1)<x_min || s_img(j,1)>x_max || s_img(j,2)<y_min || s_img(j,2)>y_max)
        continue;
    else
    [~,x_near] = min(abs(x_detect_spot(:) - s_img(j,1)));
    [~,y_near] = min(abs(y_detect_spot(:) - s_img(j,2)));
    img_distort_temp(y_near,x_near,1)  = img_distort_temp(y_near,x_near,1) + s_img(j,3);
    img_distort_temp(y_near,x_near,2)  = img_distort_temp(y_near,x_near,2) + s_img(j,4);
    img_distort_temp(y_near,x_near,3)  = img_distort_temp(y_near,x_near,3) + s_img(j,5);
    pixel_ray_count(y_near,x_near,1)  = pixel_ray_count(y_near,x_near,1) + 1;
    pixel_ray_count(y_near,x_near,2)  = pixel_ray_count(y_near,x_near,3) + 1;
    pixel_ray_count(y_near,x_near,3)  = pixel_ray_count(y_near,x_near,3) + 1;
    end
end
img_distort = img_distort_temp./pixel_ray_count;
img_distort(isnan(img_distort)) = 0;
img_distort_uint = uint8(img_distort);
figure(2);
subplot(1,1,1);imshow(img_distort_uint);title(['畸变图像（',num2str(psf_x_num),' x ',num2str(psf_y_num),' psf 网格）']);

toc;
% delete(gcp('nocreate'));

function I_color = I_add(a,x_min,x_max,x_pixels,y_min,y_max,y_pixels)
s_img = a.s_back;
nray = size(s_img,1);
x_detect_spot = linspace(x_min,x_max,x_pixels);
y_detect_spot = linspace(y_min,y_max,y_pixels);
I_color = zeros(x_pixels,y_pixels);
for j = 1:nray
    if(s_img(j,1)<x_min || s_img(j,1)>x_max || s_img(j,2)<y_min || s_img(j,2)>y_max)
        continue;
    else
    [~,x_near] = min(abs(x_detect_spot(:) - s_img(j,1)));
    [~,y_near] = min(abs(y_detect_spot(:) - s_img(j,2)));
    I_color(y_near,x_near)  = I_color(y_near,x_near) + 1;
    end
end
nray_field = sum(I_color,'all');
I_color = flipud(I_color)/nray_field;
end
