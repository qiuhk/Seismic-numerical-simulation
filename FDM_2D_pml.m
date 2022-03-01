clear
clc
tic
%%%�����Ϊ2n�ײ�֣��м�ΪPML��������NX�ż��㣬��sum_Nx = PML_n*2+n*2+Nx;ֻҪNX������Ķ��ǲ�Ҫ��
%%����pml��˼·������pml��������ӣ�Ȼ����뵽���¹������������
%% ��������
dt = 0.001; 
dh = 10; 
nx = 200;%nx�и���
n=1;   %�ռ�2n�ײ��
ny = 200;
nt = 800;
f = 10; %f<25��ɫɢ
V=4000;
slice = 20;%ÿ���100��������
slice_num = nt/slice;
%% PML��ʼ��
PML_n=20;
sum_Nx = PML_n*2+n*2+nx;
sum_Ny = PML_n*2+n*2+ny;
A1 =zeros(sum_Nx,sum_Ny);       %˥������
A2 =zeros(sum_Nx,sum_Ny);       %˥������
A3 =zeros(sum_Nx,sum_Ny);       %˥������
A4 =zeros(sum_Nx,sum_Ny);       %˥������
A =zeros(sum_Nx,sum_Ny);       %˥������

%% ��Դ
t = 1:nt;
st = 1:nt;
t0 = 100;
s_t = (1-2*(pi*f*dt*(t-t0)).^2).*exp(-(pi*f*dt*(t-t0)).^2);%Դ
%������ʼ��
P_current = zeros(sum_Nx,sum_Ny);
P_next = zeros(sum_Nx,sum_Ny);
P_past = zeros(sum_Nx,sum_Ny);
P_point = zeros(nt);
P_slice =zeros(sum_Nx,sum_Ny,slice_num);
slice_count = 1;

%% ������˥������
B = 100;  % ˥����������
for i=n+1:sum_Nx-n
    for j=n+1:sum_Ny-n  %ȥ���������
        if j<=PML_n+n
            A1(i,j) = B*(1-cos(pi*(PML_n+n+1-j)/(2*PML_n)));  %������
        elseif j>n+PML_n+ny
            A2(i,j) = B*(1-cos(pi*(j-(n+PML_n+ny))/(2*PML_n)));  %������
        end
        
        if i<=PML_n+n
            A3(i,j) = B*(1-cos(pi*(PML_n+n+1-i)/(2*PML_n)));  %������
        elseif i>n+PML_n+nx
            A4(i,j) = B*(1-cos(pi*(i-(n+PML_n+nx))/(2*PML_n)));  %������
        end
    end
end
A=A1+A2+A3+A4;%��ǰ��ö�Ӧ�����˥�����ӣ�ֱ���������¾��У�AԽ����������ԽС
clear A1 A2 A3 A4
%%%%˥������ͼ
figure(1)
contourf(A);
colorbar;
title('˥������');
%% ��ʼ����
start_time = clock;
for T= 1:nt
    %��Դ
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
    %����
    P_past = P_current;
    P_current = P_next;
  %  P_point(T)=P_current(100,300); %�첨��
    if(mod(T,slice)==0)
        run_time = etime(clock,start_time);
        fprintf('step=%d,total=%d,�ۼƺ�ʱ%.2fs\n',T,nt,run_time);
        P_slice(:,:,slice_count) = P_next;
        slice_count = slice_count + 1;
    end
end
%% ������ͼ
fmat=moviein(slice_num);
filename = 'FDM_2D_pml.gif';
figure(1)
for II = 1:slice_num
    figure(2)
    imagesc(P_slice(:,:,II));
    shading interp;
    axis tight;
    set(gca,'yDir','reverse'); %��תy��
    str_title = ['FDM 2D pml  t=',num2str(dt*II*100),'s'];
    title(str_title)
    colorbar
    caxis([-5*10^(-7) 5*10^(-7)])%��ɫ��Χ
    drawnow; %ˢ����Ļ
    F = getframe(gcf);%����ͼ����ΪӰƬ֡
    I = frame2im(F); %����ͼ������
    [I, map] = rgb2ind(I, 256); %��rgbת��������ͼ��
    if II == 1
        imwrite(I,map, filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(I,map, filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    fmat(:,II)=getframe;
end
%movie(fmat);%�ٷ�һ��
