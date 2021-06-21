% Channels positions in NeuroNexus probe 4 shanks 8 channels per shank
% dimensions in micrometres, shanks ordererd from left to right

% Basic array
array = zeros(8,2);
array(1:8,1)= zeros(1,8); 
array(1:8,2)= 50:50:400;
figure
plot(array(:,1),array(:,2),'*');

% Shank 1
Shank_1 = array;
Shank_1(:,1) = Shank_1(:,1) - 300; % move x_pos left of 300um
% Shank 2
Shank_2 = Shank_1;
Shank_2(:,1) = Shank_2(:,1) + 200; % move x_pos right of 200um
% Shank 3
Shank_3 = Shank_1;
Shank_3(:,1) = Shank_3(:,1) + 400; % move x_pos right of 200um
% Shank 4
Shank_4 = Shank_1;
Shank_4(:,1) = Shank_4(:,1) + 600; % move x_pos right of 200um

figure
plot(Shank_1(:,1),Shank_1(:,2),'*');
hold all
plot(Shank_2(:,1),Shank_2(:,2),'*');
hold all
plot(Shank_3(:,1),Shank_3(:,2),'*');
hold all
plot(Shank_4(:,1),Shank_4(:,2),'*');
