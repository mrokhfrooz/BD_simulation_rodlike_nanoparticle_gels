clc
clear
b=1000*1e-9;
dataX = importdata('save_qx0.txt') ;
dataY = importdata('save_qy0.txt') ;
dataZ = importdata('save_qz0.txt') ;


dataX= dataX*10^9;
dataY= dataY*10^9;
dataZ= dataZ*10^9;


particleid=2;% particle id to show
figure ()
dashbox2(b*10^9);

ax = gca;
ax.FontSize = 20;  % Font Size
axis equal;
for j = 1:20000   
        plot3(dataX(1,j),dataY(1,j),dataZ(1,j),'ko')  
end



