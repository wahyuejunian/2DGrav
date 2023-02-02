% By: Wahyu Eko Junian
% 2D Gravity inversion  (Talwani)

clc; close all; clear

%__________________________________________%
G = 6.673*10^-11; %konstanta gravitasi umum
z0=0; %topografi dianggap konstan (0)
%__________________________________________%
dx=20; % besar grid arah x
dz=40; % besar grid arah z
x=0:50:2000; % panjang x-position
z=0:50:1000; % kedalaman z-position
node = 4; % banyak titik sudut polygon
%__________________________________________%
xc=0:dx:2000; % banyak cell arah x
zc=0.1:dz:1000; % banyak cell arah z
%__________________________________________%
density=zeros(length(zc)-1,length(xc)-1);
density(8:13,60:80)=1000; % Benda horizontal
density(4:16,17:22)=1000; % Benda vertical
%__________________________________________%
gobs=zeros(1,length(x));
kernel=zeros(length(x),(length(xc)-1)*(length(zc)-1));
for ii = 1:length(x)
    k=1;
    sum=0;
    for i = 1:length(xc)-1
        for j = 1:length(zc)-1
            % titik sudut cell pada xc dan zc
            xm = [xc(i);xc(i+1);xc(i+1);xc(i)];
            zm = [zc(j);zc(j);zc(j+1);zc(j+1)];
            
            sum=sum+2*G*density(j,i)*Talwani(x(ii),z0,xm,zm,node);% respons forward
            kernel(ii,k)=2*G*Talwani(x(ii),z0,xm,zm,node); % matrix kernel
            k=k+1;
        end
    end
    gobs(ii) = sum; % respons
end
%__________________________________________%
% proses inversi
model = kernel'*inv(kernel*kernel')*gobs'; % m=G'*inv(GG')*d
%__________________________________________%
gcal=kernel*model; %respons inversi d=G.m

model_pred=reshape(model,size(density));

%% visualisasi
fh = figure();
subplot(3,1,1) % gambar respons forward vs inversi
hold on
plot (x,gobs,'or','Markerfacecolor','r');
plot(x,gcal,'-b','LineWidth',2);
legend ('Synthetic Data', 'Predicted Data')
title ('Gravity Anomaly')
ylabel('\bf mgal');

subplot(3,1,2) % gambar model inversi
colormap jet
imagesc(xc(1:100),zc(1:24),model_pred)
ylabel('\bf Depth (m)');
set(gca,'ydir','reverse');
title('Invertion Model (unconstrained strategy)')


subplot(3,1,3) % gambar model synthetic
% colormap jet
for i=1:length(xc)-1
    for j= 1:length(zc)-1
        xfill = [xc(i);xc(i+1);xc(i+1);xc(i)];
        zfill = [zc(j);zc(j);zc(j+1);zc(j+1)];
        val = density(j,i);
        patch(xfill,zfill,val);  
    end
end
xlim([0 2000])
ylim([0 max(zfill)])
title('Synthetic Model (density 1 g/cc)')
xlabel('\bf Distance (m)');
ylabel('\bf Depth (m)');
set(gca,'ydir','reverse');


fh.WindowState = 'maximized';
print('-dpng','GravInv_Matlab','-r300');


