clear
clc
tic
%% ��������
dt = 0.001; 
dh = 10; 
nx = 400;%nx����
ny = 400;
nn=100;  %���ղ�
nt =900;
f = 10; %f<25��ɫɢ
V=4000;
slice = 20;%ÿ���100��������
slice_num = nt/slice;
%% ��ʼ��
%% ��Դ
t = 1:nt;
st = 1:nt;
t0 = 100;
s_t = (1-2*(pi*f*dt*(t-t0)).^2).*exp(-(pi*f*dt*(t-t0)).^2);%Դ
%������ʼ��
P_current = zeros(nx,ny);
P_next = zeros(nx,ny);
P_past = zeros(nx,ny);
ddp=zeros(nx,ny);
xi=ones(nx,ny);%���ղ�
P_point = zeros(nt);
P_slice =zeros(nx,ny,slice_num);
slice_count = 1;

%% ��ʼ����
start_time = clock;
for T= 1:nt
    %��Դ
    P_current(nx/2,ny/2) = P_current(nx/2,ny/2) + dt^2*s_t(T);
    for j=2:nx-1
        for k=2:ny-1
           ddp(j,k)=(P_current(j+1,k)-2*P_current(j,k)+P_current(j-1,k))/dh/dh...
               +(P_current(j,k+1)-2*P_current(j,k)+P_current(j,k-1))/dh/dh;
        end
    end
    P_next=V*V*dt*dt*ddp+2*P_current-P_past;
    %���ղ� Ч������
    for j=1:nx
        for k=1:nn
           xi(j,k)=exp(-1*(0.0005*(nn-k))^2);%��
        end
    end
    for j=1:nx
        for k=ny-nn:ny
           xi(j,k)=exp(-1*(0.0005*(ny-nn-k))^2);
        end
    end
    for k=1:ny  %kΪ��
        for j=1:nn
           xi(j,k)=exp(-1*(0.0005*(nn-j))^2);%��
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
    P_point(T)=P_current(100,300); %�첨��
    if(mod(T,slice)==0)
        run_time = etime(clock,start_time);
        fprintf('step=%d,total=%d,�ۼƺ�ʱ%.2fs\n',T,nt,run_time);
        P_slice(:,:,slice_count) = P_next;
        slice_count = slice_count + 1;
    end
end
%% ������ͼ
fmat=moviein(slice_num);
filename = 'FDM_2D.gif';
figure(1)
for II = 1:slice_num
    imagesc(P_slice(:,:,II));
    shading interp;
    axis tight;
    grid on;
    set(gca,'yDir','reverse'); %�̶�colorbar
    str_title = ['FDM 2D  t=',num2str(dt*II*100),'s'];
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
