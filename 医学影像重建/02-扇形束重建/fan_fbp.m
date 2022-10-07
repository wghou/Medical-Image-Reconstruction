clear;
clc;

I = phantom(256);   % 获得脑部影像
% I = imread("panda.png");
% I = rgb2gray(I);
beta = 0:1:359;     % 探测器的旋转角
 

% === step1:获得扇束投影 === %
[m,n] = size(I);
% bug在这！！ 2*m -> 20*m
d = 20 * m;    % 射线源与旋转中心之间的距离
[F, detector, angle] = fanbeam(I, d, 'FanSensorGeometry', 'line');



% === step2:滤波前加权weight1 === %
detector_num = length(detector);    % 记录探测器的数量    

tensor_d = ones(detector_num,1).*d;   % 将d扩展为列向量，方便运算
weight1 = tensor_d ./ (sqrt(tensor_d.*tensor_d + detector.*detector));% 计算权重weight1

% 对投影加权
for i = 1:length(beta)
    F(:,i) = F(:,i) .* weight1;  
end


tic;

% === step3:滤波 === %
% % >> 方法一: 时域卷积
% % 构造滤波器filter
% t = linspace(-n/2, n/2-1, n);
% filter = 0.0085 * (sinc(t)/2-sinc(t/2).^2/4);
% % 调整投影位置,首位垫m行
% R = F;
% F = zeros(detector_num+2*m, 360);
% F(m+1:m+detector_num, :) = R;
% for N = 1:360
%     temp(:, N) = conv(F(:,N), filter);
% end
% % 选取中间的有效部分作为卷积结果(m行之前的留白+m/2行的脉冲函数)
% F = temp(3*m/2:3*m/2+detector_num, :);

% >> 方法二: 频域乘积
% 得到投影F的傅立叶变换f
f = fft(F, detector_num);
% 这里采用Ram-Lak滤波器
filter = 2 * [0:round(detector_num/2-1), detector_num/2:-1:1]' / detector_num;
f_filtered = f .* filter;     % 滤波：用投影结果点乘filter
% 傅里叶逆变换
F = ifft(f_filtered);
F = real(F);  % 取实部

toc;

% % 使用ifanbeam检查一下效果, F正常
% I1 = ifanbeam(F,d, 'FanSensorGeometry', 'line');
% imshow(I1, []);



% 检查: 错误与探测器数量之间的关系
% === step4:加权反投影 === %
I1 = zeros([m,n]);

for b = beta
    r = deg2rad(b); % 角度转为弧度,便于使用sin/cos
    % 设当前待重建点坐标为(x,y)
    for y = 1:m
        for x = 1:n
            % 计算s': 相似三角形
            s1 = d * ((x-n/2)*cos(r)+(y-m/2)*sin(r)) / (d+(x-n/2)*sin(r)-(y-m/2)*cos(r));
            % 确保s'未超出探测器范围
            if (s1>=-(detector_num-1)/2) && (s1<=(detector_num-1)/2)
                % 取整,得到该点(亦即s')对应点detector序号
                s1_int = round(s1);     
                % 重建该点
                point = F(s1_int+(detector_num-1)/2+1, b+1);
                % 积分反投影 (离散积分用求和表示): 注意matlab中纵轴的方向
                I1(m-y+1,x) = I1(m-y+1,x) + 1/(2*power(U(d,x-m/2,y-m/2,r),2)) * point;
            end
        end
    end
end    


imshow(I1,[]);




% === func:计算公式中的 u === %
function u = U(d,x,y,beta)
    u = (d+x*sin(beta)-y*cos(beta)) / d;
end







