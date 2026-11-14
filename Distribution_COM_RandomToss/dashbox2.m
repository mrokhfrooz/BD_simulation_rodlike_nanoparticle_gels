function x =plotBox2(b)

x = [0 b b 0 ] ;%0 10 10 breadth-10 breadth-10 10];
y = [0 b];
z = [0 0 b b ] ;%0 10 10 breadth-10 breadth-10 10];
m = length(x);
xr=[x x(1)];
yr1=y(1)*ones(m+1,1);
zr=[z z(1)];
plot3(xr,yr1,zr,'--k');
hold on
yr2=y(2)*ones(m+1,1);
plot3(xr,yr2,zr,'--k');
hold on
for i=1:1:m
plot3([x(i) x(i)],y,[z(i) z(i)],'--k');    
end

xlabel('X (nm)');
ylabel('Y (nm)');
zlabel('Z (nm)');
view(-37.5,30)
%view (2); % https://www.mathworks.com/help/matlab/ref/view.html
%axis equal;

hold on
