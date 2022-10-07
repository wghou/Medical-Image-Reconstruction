clc;
clear;
close all;

% ==== step1:计算每个投影g(p,theta)的一维傅里叶变换G(w,theta) ==== %
I = imread("panda.png");
I = rgb2gray(I);    % 转为灰度图片

theta = 0:179;
g = radon(I,theta);   % 做0~179度共计180个角度的投影，得到g(p,theta)
width = length(g); 

G = fft(g, width);   % 一维快速傅里叶变换，得到G(w,theta)。注意：该结果为复数


% ==== step2:将G(w,theta)乘以斜坡滤波器|w| ==== %
% 这里采用Ram-Lak滤波器
filter = 2*[0:round(width/2-1), width/2:-1:1]' / width;
G_filtered = G .* filter;     % 滤波：用投影结果点乘filter


% ==== step3:对step2的结果做一维傅里叶逆变换 ==== %
g_filtered = ifft(G_filtered);
g_filtered = real(g_filtered);  % 取实部


% ==== step4:将所有处理后的结果叠加，最终得到原图像。 ==== %
result = iradon(g_filtered, theta, 'none'); % 不使用iradon自带的滤波器，直接反投影
imshow(result, []), title("FBP result");


% 参照1:直接反投影，无滤波
figure; imshow(iradon(g, theta, 'none'), []), title("参照1: 直接反投影，无滤波");
% 参照2:iradon自带"Ram-Lak"滤波
figure; imshow(iradon(g, theta, 'Ram-Lak'), []), title("参照2: iradon自带Ram-Lak滤波"); 

