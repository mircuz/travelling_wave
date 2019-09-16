close all
%clear all

%results = importfile("results.txt");
%fprintf("Data read\n");

nx= 100;
ny= 100;
res= reshape(results',nx,ny,[]);

fprintf("Data reshaped, ready for plot\n");

x = linspace(0,1,nx);
y = linspace(0,1,ny);

for k = 1:size(res,3)
    surf(x,y,res(:,:,k));
    zlim([-1 1]);
    colormap(copper);
    %view(2);
    zlim([-1 1]);
    drawnow
    pause(0.01);
end