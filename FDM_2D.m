clear
clc
tic
%% 参数设置
dt = 0.001; 
dh = 10; 
nx = 400;%nx个点
ny = 400;
nn=100;  %吸收层
nt =900;
f = 10; %f<25无色散
V=4000;
slice = 20;%每间隔100步保存结果
slice_num = nt/slice;
%% 初始化
%% 震源
t = 1:nt;
st = 1:nt;
t0 = 100;
s_t = (1-2*(pi*f*dt*(t-t0)).^2).*exp(-(pi*f*dt*(t-t0)).^2);%源
%变量初始化
P_current = zeros(nx,ny);
P_next = zeros(nx,ny);
P_past = zeros(nx,ny);
ddp=zeros(nx,ny);
xi=ones(nx,ny);%吸收层
P_point = zeros(nt);
P_slice =zeros(nx,ny,slice_num);
slice_count = 1;

%% 开始计算
start_time = clock;
for T= 1:nt
    %加源
    P_current(nx/2,ny/2) = P_current(nx/2,ny/2) + dt^2*s_t(T);
    for j=2:nx-1
        for k=2:ny-1
           ddp(j,k)=(P_current(j+1,k)-2*P_current(j,k)+P_current(j-1,k))/dh/dh...
               +(P_current(j,k+1)-2*P_current(j,k)+P_current(j,k-1))/dh/dh;
        end
    end
    P_next=V*V*dt*dt*ddp+2*P_current-P_past;
    %吸收层 效果不好
    for j=1:nx
        for k=1:nn
           xi(j,k)=exp(-1*(0.0005*(nn-k))^2);%左
        end
    end
    for j=1:nx
        for k=ny-nn:ny
           xi(j,k)=exp(-1*(0.0005*(ny-nn-k))^2);
        end
    end
    for k=1:ny  %k为列
        for j=1:nn
           xi(j,k)=exp(-1*(0.0005*(nn-j))^2);%上
        end
    end
     for k=1:ny
        for j=ny-nn:ny
           xi(j,k)=exp(-1*(0.0005*(nx-nn-j))^2);
        end
     end
    P_next=P_next.*xi;
    
    P_past = P_current;
    P_current = P_next;
    P_point(T)=P_current(100,300); %检波器
    if(mod(T,slice)==0)
        run_time = etime(clock,start_time);
        fprintf('step=%d,total=%d,累计耗时%.2fs\n',T,nt,run_time);
        P_slice(:,:,slice_count) = P_next;
        slice_count = slice_count + 1;
    end
end
%% 制作动图
fmat=moviein(slice_num);
filename = 'FDM_2D.gif';
figure(1)
for II = 1:slice_num
    imagesc(P_slice(:,:,II));
    shading interp;
    axis tight;
    grid on;
    set(gca,'yDir','reverse'); %固定colorbar
    str_title = ['FDM 2D  t=',num2str(dt*II*100),'s'];
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
