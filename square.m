clc;clear;close all;
%%
%��ȡ���ݣ���ͼ�������ٶȡ�ƽ���ٶȺͼ��ٶ�
Data = importdata("E:\OneDrive - mail.ustc.edu.cn\LABVIEW\ʵ������0112\02\0201.csv");
data_1 = Data(:,1:2);
[M,N] = size(data_1);

x = [-0.035,0.035,0.035,-0.035,-0.035];
y = [-0.035,-0.035,0.035,0.035,-0.035];
figure;
plot(Data(:,1),Data(:,2),'b');
axis equal;
ylim([-0.055 0.055]);
xlim([-0.055 0.055]);
hold on;
plot(x,y,'r');
axis equal;
axis off;
print(gcf, 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\ʵ������0112\02\0201.bmp','-r600','-dbmp');
figure;plot(x,y,'r');
axis equal;
ylim([-0.055 0.055]);
xlim([-0.055 0.055]);
axis off;
print(gcf, 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\ʵ������0112\02\0.bmp','-r600','-dbmp');

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
%�ٶȷ�ֵ�㼰��ͼ
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
print(gcf, 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\ʵ������0112\02\v_0201.bmp','-r600','-dbmp')
%%
%ƽ���켣ƫ��ֵ���غ϶�
num_similar = 0;
distance = [];
for i = 1:M
    angle_1 = atan2(data_1(i,2),data_1(i,1));
    angle = rad2deg(angle_1);
    if(angle < 45 && angle > -45)
        distance(i) = abs(data_1(i,1) - 0.035);
    end
    if(angle < 135 && angle > 45)
        distance(i) = abs(data_1(i,2) - 0.035);
    end
    if(angle <= 180 && angle > 135 || angle < -135 && angle >= -180)
        distance(i) = abs(data_1(i,1) + 0.035);
    end
    if(angle < -45 && angle > -135)
        distance(i) = abs(data_1(i,2) + 0.035);
    end
    if(angle == 45 || angle == 135 || angle == -45 || angle == -135)
        dx = abs(data_1(i,1)) - 0.035;
        dy = abs(data_1(i,2)) - 0.035;
        distance(i) = sqrt(dx*dx + dy*dy);
    end
        if(distance(i) <= 0.001)
    num_similar = num_similar + 1;
    end
end
coindence = num_similar/M;
sumofd = sum(distance);
d = sumofd / M;

%%
I  = imread("E:\OneDrive - mail.ustc.edu.cn\LABVIEW\ʵ������0112\02\0201.bmp");
Ic = imread("E:\OneDrive - mail.ustc.edu.cn\LABVIEW\ʵ������0112\02\0.bmp");
I_1  = rgb2gray(I);
Ic_1 = rgb2gray(Ic);
[M,N] = size(I_1);
%%
R=I(:,:,1); %��ȡ��ɫ����
G=I(:,:,2); %��ȡ��ɫ����
B=I(:,:,3); %��ȡ��ɫ����
I_out = I;

I_2 = zeros(M,N);
Ic_2 = zeros(M,N);
for i = 1:M
    for j = 1:N
        if(I_1(i,j) == 255)
            I_2(i,j) = 1;
            Ic_2(i,j) = 1;
        end
    end
end
L = bwlabel(I_2,8);  %��ȡ��ͨ��
Lc = bwlabel(Ic_2,8);

area = 0;
squareArea = 0;
for i = 1:M
    for j = 1:N
        if(L(i,j) ~= 1 && L(i,j) ~= 3  && L(i,j) ~= 0)
            area  = area + 1;
            R(i,j) = 255;
            G(i,j) = 255;
            B(i,j) = 0;
        end
        if(Lc(i,j) == 3 )
            squareArea = squareArea + 1;
        end
    end
end
I_out(:,:,1) = R;
I_out(:,:,2) = G;
I_out(:,:,3) = B;
S = 0.07*0.07*area/squareArea;
figure;imshow(I_out);
imwrite(I_out,'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\ʵ������0112\02\0201_out.bmp');
%%
%������
fprintf('��ֵ��������%d\n',numofpks);
fprintf('ƽ���ٶȣ�%d\n',V_avg);
fprintf('ƽ�����ٶȣ�%d\n',A_avg);
fprintf('ƽ���켣ƫ�%d\n',d);
fprintf('�켣�غ϶ȣ�%f\n',coindence);
fprintf("�����ι켣�ཻ�γɵ����Ϊ��%f\n",S);