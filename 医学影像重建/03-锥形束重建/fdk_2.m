close all;
clear;
clc;


% config
n = 128;
beta = deg2rad(0:359);
beta_num = length(beta);
d = 10000; % 源到旋转中心


% === step1: 投影 === %
% 导入shepp-logan-3d和其投影
load("shepp_logan_3d_128.mat");     % I
load("proj_shepp_logan_3d_128.mat");    % P


% === step2: 加权 === %
detector = -n/2+1 : n/2;   % 探测器的实际位置
detector = detector';
% xi = ξ, 锥束投影中距中平面的距离
for xi = detector
    tensor_d = ones(n,1) .* d;
    tensor_xi = ones(n,1) .* xi;
    weight = tensor_d ./ (sqrt(tensor_d.*tensor_d + detector.*detector + tensor_xi.*tensor_xi));

    % 对投影加权
    for i = 1:beta_num
        P(xi+n/2,:,i) = P(xi+n/2,:,i) .* weight;
    end
end



% === step3: 滤波 === %
% <<<< 方法一: 频域滤波 >>>>
temp = zeros(n, n, beta_num);  % 用于暂存滤波结果
% 这里采用Ram-Lak滤波器
filter =  2 * [1:n/2, n/2:-1:1]' / n;
% 滤波: 频域乘积|w|
for z = 1:n
    f = fft(P(z,:,:),n);
    f = squeeze(f);
    f_filtered = f .* filter;
    temp(z,:,:) = ifft(f_filtered);
end
P = real(temp);  % 取实部

% % <<<< 方法二: 时域卷积 >>>>
% % 构造滤波器filter
% t = linspace(-n/2+1, n/2, n);
% filter = 5 * (sinc(t)/2-sinc(t/2).^2/4);
% % 调整投影位置,首位垫n行
% R = P;
% P = zeros(n, n+2*n, 360);
% P(:,n+1:n+n,:) = R;
% for b = 1:360
%     for z = 1:n
%         temp(z,:,b) = conv(P(z,:,b),filter);
%     end
% end
% % 选取中间的有效部分作为卷积结果(n行之前的留白+n/2行的脉冲函数)
% P = temp(:,3*n/2:3*n/2+n,:);

% % 检查滤波效果, 正常
% for i = 1:10:360
%     figure;
%     imshow(P(:,:,i),[]);
% end



% === step4: 反投影 === %
I = zeros([n,n,n]);
for b = 0:359
    r = deg2rad(b); % 角度转为弧度,便于使用sin/cos
    % 设当前待重建点坐标为(x,y,z)
    for y = 1:n     
        for x = 1:n
            for z = 1:n     % 对于任意ξ, 系数d^2/(d-s)^2都相同
                s = -(x-n/2)*sin(r) + (y-n/2)*cos(r);
                t1 = (x-n/2) *cos(r) + (y-n/2) *sin(r); % t_tilde = t, 即原来的s'
                z1 = z * d/(d-s);
                % 确保t1未超出探测器范围
                if (t1>-n/2) && (t1<=n/2) && (z1>=1) && (z1<=n)
                    % 取整,得到该点对应点detector序号
                    t1_int = ceil(t1);  
                    z1_int = ceil(z1);
                    % 重建该点
                    point = P(z1_int, t1_int+n/2, b+1);

                    % 积分反投影 (离散积分用求和表示): 注意matlab中纵轴的方向
                    I(n-x+1, y, z) = I(n-x+1, y, z) + d*d / (2*power((d-s),2)) * point;
                end
            end
        end
    end
end


for i = 1:10
% for i = 1:4:128
    subplot(2,5,i);
    imshow(I(:,:,i+n/2-5),[]);
end