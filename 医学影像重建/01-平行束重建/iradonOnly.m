clear;
close all;
% 读取图片，并将图片改为单通道
pic = imread("panda.png"); 
pic = rgb2gray(pic);  
% 投影——radon变换
theta = 0:179;
proj = radon(pic, theta);

figure;
imshow(iradon(proj, theta, 'none'), []);
figure;
imshow(iradon(proj, theta),[]);


