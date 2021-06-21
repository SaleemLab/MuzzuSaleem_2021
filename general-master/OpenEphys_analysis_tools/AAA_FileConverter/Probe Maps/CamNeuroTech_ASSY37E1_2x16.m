% Channels positions in NeuroNexus probe 4 shanks 8 channels per shank
% dimensions in micrometres, shanks ordererd from left to right

% Basic array
array = zeros(16,2);
array(1:16,1)= zeros(1,16); 
array(1:16,2)= [0 linspace(25,300,15)];
figure
plot(array(:,1),array(:,2),'*');

electrode_dim = [11,15]; % base by height
% Shank 1
Shank_1 = array;
g = 1;
for i=2:2:16
    Shank_1(i,1)= -(g)*electrode_dim(1)/4;
    g=g+1;
end
g = 2;
for i=3:2:15
    Shank_1(i,1)= (g)*electrode_dim(1)/4;
    g=g+1;
end
Shank_1(:,1) = Shank_1(:,1) - 125; % move x_pos left of 125um
% Shank 2
Shank_2 = Shank_1;
% suggested fix for problem with 2nd shank (Jonas)
% Shank_2(:,1) = 0;
g = 1;
for i=2:2:16
    Shank_2(i,1)= -(g)*electrode_dim(1)/4;
    g=g+1;
end
g = 2;
for i=3:2:15
    Shank_2(i,1)= (g)*electrode_dim(1)/4;
    g=g+1;
end
Shank_2(:,1) = Shank_2(:,1) + 125; % move x_pos right of 125um
Shank_2(1,1) = Shank_2(1,1) + 125; % move x_pos right of 125um

figure
plot(Shank_1(:,1),Shank_1(:,2),'*');
hold all
plot(Shank_2(:,1),Shank_2(:,2),'*');
