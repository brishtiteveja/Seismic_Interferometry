% code to generate station.csv to use in this matlab file
% st_sum_matrix = matrix(unlist(st_sum), ncol = 20, byrow = FALSE)
% write.table(st_sum_matrix, file="./station.csv", sep=",", col.names = F, row.names = F)

addpath(genpath('./crewes'))
clc ; clear ; close all;

M = dlmread('station.csv',',',1,0);

n = 32000
m = n / 4
mm = m/2 - 1

M = M(1:mm,:) ;

St1 = M(:,1) ;

St1 = St1 .*tukeywin(length(St1),10);

dt = 0.004;
t = (0 : length(St1)-1).*dt;

% number of station
num_station = 20

dist_station = 100

% taper perc
taper_perc = 5

% frequency filters
f_min = 1
f_max = 8

% plot(St1)
% 
% [spec,f]= fftrl(St1,t);
% 
% figure
% plot(f,abs(spec))
% 
% 
% trout=filtf(St1,t,15,20);
% 
% figure
% plot(trout)
% 
% [spec,f]= fftrl(trout,t);
% 
% figure
% plot(f,abs(spec))
[m n] = size(M);

for i = 1 : n
    M(:,i) =  M(:,i) .*tukeywin(m,taper_perc);
end

for i = 1 : n
    M(:,i) = filtf( M(:,i),t,f_min,f_max);
end

x = (0:num_station - 1).*dist_station ;

close all
%imagesc(x,t,M)

[dout] = AGCgain(M,dt,2.5,1);


[X,T] = meshgrid(x,t) ;

nx = 0:0.5:((num_station - 1)*100);
[Xx,NT] = meshgrid(nx,t) ;

new_res = interp2(X,T,dout,Xx,NT);
h = figure

imagesc(nx,t,new_res)

print(h,'structure.pdf','-dpdf')
