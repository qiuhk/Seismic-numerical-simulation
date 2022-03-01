clear
clc
tic
%%%最外层为2n阶差分，中间为PML，最里面NX才计算，共sum_Nx = PML_n*2+n*2+Nx;只要NX，其余的都是不要的
%%设置pml的思路是先求pml区域的因子，然后带入到更新过程中运算就行
%% 参数设置
dt = 0.001; 
dh = 10; 
nx = 200;%nx行个点
n=1;   %空间2n阶差分
ny = 200;
nt = 800;
f = 10; %f<25无色散
V=4000;
slice = 20;%每间隔100步保存结果
slice_num = nt/slice;
%% PML初始化
PML_n=20;
sum_Nx = PML_n*2+n*2+nx;
sum_Ny = PML_n*2+n*2+ny;
A1 =zeros(sum_Nx,sum_Ny);       %衰减因子
A2 =zeros(sum_Nx,sum_Ny);       %衰减因子
A3 =zeros(sum_Nx,sum_Ny);       %衰减因子
A4 =zeros(sum_Nx,sum_Ny);       %衰减因子
A =zeros(sum_Nx,sum_Ny);       %衰减因子

%% 震源
t = 1:nt;
st = 1:nt;
t0 = 100;
s_t = (1-2*(pi*f*dt*(t-t0)).^2).*exp(-(pi*f*dt*(t-t0)).^2);%源
%变量初始化
P_current = zeros(sum_Nx,sum_Ny);
P_next = zeros(sum_Nx,sum_Ny);
P_past = zeros(sum_Nx,sum_Ny);
P_point = zeros(nt);
P_slice =zeros(sum_Nx,sum_Ny,slice_num);
slice_count = 1;

%% 余弦型衰减因子
B = 100;  % 衰减幅度因子
for i=n+1:sum_Nx-n
    for j=n+1:sum_Ny-n  %去除差分区域
        if j<=PML_n+n
            A1(i,j) = B*(1-cos(pi*(PML_n+n+1-j)/(2*PML_n)));  %左区域
        elseif j>n+PML_n+ny
            A2(i,j) = B*(1-cos(pi*(j-(n+PML_n+ny))/(2*PML_n)));  %右区域
        end
        
        if i<=PML_n+n
            A3(i,j) = B*(1-cos(pi*(PML_n+n+1-i)/(2*PML_n)));  %下区域
        elseif i>n+PML_n+nx
            A4(i,j) = B*(1-cos(pi*(i-(n+PML_n+nx))/(2*PML_n)));  %上区域
        end
    end
end
A=A1+A2+A3+A4;%提前求好对应区域的衰减因子，直接拿来更新就行，A越是里面区域越小
clear A1 A2 A3 A4
%%%%衰减因子图
figure(1)
contourf(A);
colorbar;
title('衰减因子');
%% 开始计算
start_time = clock;
for T= 1:nt
    %加源
    P_current(sum_Nx/2,sum_Ny/2) = P_current(sum_Nx/2,sum_Ny/2) + dt^2*s_t(T);
    for i=n+1: sum_Nx-n
        for j=n+1: sum_Ny-n
           pml1 = (2-A(i,j)*A(i,j)*dt*dt)/(1+A(i,j)*dt)*P_current(i,j);
           pml2 = (1-A(i,j)*dt)/(1+A(i,j)*dt)*P_past(i,j);
           ddp=(P_current(i+1,j)-2*P_current(i,j)+P_current(i-1,j))/dh/dh...
               +(P_current(i,j+1)-2*P_current(i,j)+P_current(i,j-1))/dh/dh;
           P_next(i,j) = pml1 - pml2 + (V*V*dt*dt)/(1+A(i,j)*dt)*ddp;
        end
    end
    %更新
    P_past = P_current;
    P_current = P_next;
  %  P_point(T)=P_current(100,300); %检波器
    if(mod(T,slice)==0)
        run_time = etime(clock,start_time);
        fprintf('step=%d,total=%d,累计耗时%.2fs\n',T,nt,run_time);
        P_slice(:,:,slice_count) = P_next;
        slice_count = slice_count + 1;
    end
end
%% 制作动图
fmat=moviein(slice_num);
filename = 'FDM_2D_pml.gif';
figure(1)
for II = 1:slice_num
    figure(2)
    imagesc(P_slice(:,:,II));
    shading interp;
    axis tight;
    set(gca,'yDir','reverse'); %翻转y轴
    str_title = ['FDM 2D pml  t=',num2str(dt*II*100),'s'];
    title(str_title)
    colorbar
    caxis([-5*10^(-7) 5*10^(-7)])%颜色范围
    drawnow; %刷新屏幕
    F = getframe(gcf);%捕获图窗作为影片帧
    I = frame2im(F); %返回图像数据
    [I, map] = rgb2ind(I, 256); %将rgb转换成索引图像
    if II == 1
        imwrite(I,map, filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(I,map, filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    fmat(:,II)=getframe;
end
%movie(fmat);%再放一遍
