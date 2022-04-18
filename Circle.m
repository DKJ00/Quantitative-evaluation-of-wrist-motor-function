clc;clear;close all;
%%
%读取数据，绘图，计算速度、平均速度和加速度
Data = importdata("E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0112\02\0021.csv");
data_1 = Data(:,1:2);
[M,N] = size(data_1);

r=0.045; theta=0:pi/128:2*pi;
x=r*cos(theta); y=r*sin(theta);
figure;
plot(Data(:,1),Data(:,2),'b');
ylim([-0.06 0.06]);
xlim([-0.06 0.06]);
hold on;
plot(x,y,'r');
axis equal;
axis off;
print(gcf, 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0112\02\0021.bmp','-r600','-dbmp');
figure;plot(x,y,'r');
axis equal;
ylim([-0.06 0.06]);
xlim([-0.06 0.06]);
axis off;
print(gcf, 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0112\02\0.bmp','-r600','-dbmp');
fs = 100;
t = 1/fs;
for i = 1:M-1
    dx = abs(data_1(i+1,1) - data_1(i,1));
    dy = abs(data_1(i+1,2) - data_1(i,2));
    velo(i) = sqrt(power(dx,2) + power(dy,2))/t;
end
V_avg = mean(velo);
for i = 1:M-2
    acce(i) = (velo(i+1) - velo(i))/t;
end
A_avg = mean(abs(acce(2:end)));

%%
%速度峰值点及绘图
[pks,locs] = findpeaks(velo,'minpeakdistance',1);
numofpks = 0;
for i = 1:length(pks)
    if(pks(i)>= 1.5*V_avg)
        numofpks = numofpks + 1;
    end
end
figure;
plot(velo);
hold on;
for i = 1:length(pks)
    text(locs(i),pks(i),'*','color','g');
end
hold on;
line_E = 1.5*V_avg*ones(1,length(velo));
plot(line_E,'r');
print(gcf, 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0112\02\v_0021.bmp','-r600','-dbmp')
%%
%轨迹偏差值及重合度
num_similar = 0;
distance = [];
for i = 1:M
    dx = abs(data_1(i,1));
    dy = abs(data_1(i,2));
    distance(i) = abs(sqrt(power(dx,2) + power(dy,2)) - r);
    if(distance(i) <= 0.001)
        num_similar = num_similar + 1;
    end
end
coindence = num_similar/M;
sumofd = sum(distance);
d = sumofd / M;

%%
I  = imread("E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0112\02\0021.bmp");
Ic = imread("E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0112\02\0.bmp");
I_1  = rgb2gray(I);
Ic_1 = rgb2gray(Ic);
[M,N] = size(I_1);
%%
R=I(:,:,1); %获取红色分量
G=I(:,:,2); %获取绿色分量
B=I(:,:,3); %获取蓝色分量
I_out = I;

I_2  = zeros(M,N);
Ic_2 = zeros(M,N);
for i = 1:M
    for j = 1:N
        if(I_1(i,j) == 255)
            I_2(i,j)  = 1;
            Ic_2(i,j) = 1;
        end
    end
end
L  = bwlabel(I_2,8);  %获取连通域
Lc = bwlabel(Ic_2,8);

area = 0;
circleArea = 0;
for i = 1:M
    for j = 1:N
        if(L(i,j) ~= 1 && L(i,j) ~= 6 && L(i,j) ~= 0)
            area  = area + 1;
            R(i,j) = 255;
            G(i,j) = 255;
            B(i,j) = 0;
        end
        if(Lc(i,j) ==  6  )
            circleArea = circleArea + 1;
        end
    end
end
I_out(:,:,1) = R;
I_out(:,:,2) = G;
I_out(:,:,3) = B;
S = pi*0.045*0.045*area/circleArea;
figure;imshow(I_out);
imwrite(I_out,'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0112\02\0021_out.bmp');

%%
%结果输出
fprintf('峰值点数量：%d\n',numofpks);
fprintf('平均速度：%d\n',V_avg);
fprintf('平均加速度：%d\n',A_avg);
fprintf('平均轨迹偏差：%d\n',d);
fprintf('轨迹重合度：%f\n',coindence);
fprintf("圆形轨迹相交形成的面积为：%f\n",S);