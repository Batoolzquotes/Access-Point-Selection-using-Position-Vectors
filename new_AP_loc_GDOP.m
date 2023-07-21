
% value_calc_GDOP_0601

clc
clear 
warning off

AP_num = 1;       %using 1 additional AP
AP_init_num = 3;  %3 already deployed APs
L=8;     %gray coding bit length
N=60;    %60 Reference points
M=98;    % GA chromosome length
T=300;   %GA evolution time
Pc=0.9;  %GA crossover probability        
Pm=0.02; %GA mutation probability
PT = 256;  %GA population size    %% x,z genes, set of 2 genes for each location, 256 sets in population
m = 16;    % 16x32 grid
n =32;     % = 512 grid points
AP_total = AP_num + AP_init_num; 
rand_value = normrnd(0,0,AP_total,m*n);
grid_size =1; 
EU_thr = 2; 

for i=1:1:N  %for 60 reference locations
    x1(:,i)=rand(1,AP_num)*m;  %generate 'm' random numbers with maximum 1
    x2(:,i)=uint8(x1(:,i)/m*PT); %256 RSSI readings at each RP   %%convert number to 8 bit unsigned integer, averaged over m, 256 population size
    z1(:,i) = rand(1,AP_num)*n; %generate 'n' random numbers with maximum 1 
    z2(:,i)=uint8(z1(:,i)/n*PT); %256 RSSI readings at each RP   %%convert number to 8 bit unsigned integer, averaged over n, 256 population size
    for j = 1:AP_num 
    grayCode_x(j,i,:)=num2gray(x2(j,i),L);  %convert 8-bit integers to gray coded (consecutive numbers differ by 1 bit) using function 'num2gray' defined in another m-file
    grayCode_z(j,i,:)=num2gray(z2(j,i),L);  %similarly for z-coordinate 
    end  
end
%%%%%%%%%

G_E(1) = 0;  % set GA evolution time to zero
same_num = 0; 
for t=1:1:T   
    t   %display t=1:300 with increment of 1

power_trans = 24;   %transmit power specified
nh = 3.1;   %attenuation factor of indoor environment
M_f = 2;   
loss_a = 3.1;  
N_f = 2;
loss_b = 3.1;
for i =1:m     %for 16 values of i (iterate for x-plane)
    for j = 1:n    %for 32 values of j (iterate for z-plane) 
        k = (j-1)*m +i;  
        a(k,1) = 0.5+(i-1);   % Specifying 512 Grid points (mid points of grid)
        a(k,2) = 0.5 + grid_size*(j-1);   
    end
end
RP_loc = a;   % map RP locations to grid points
y1 = value_calc_GDOP_0601(x1,z1,N,AP_num,RP_loc,AP_init_num);  %calculate GDOP from function
for i=1:1:M/2  %
    [a,b]=min(y1) ;   %find location where GDOP is minumum
    for k = 1:AP_num  
    grayCodeNew_x(k,i,:)=grayCode_x(k,b,:);  %gray code the new location
    grayCodeNew_x(k,i+M/2,:)=grayCode_x(k,b,:); %map GA string to new location
    
    
    grayCodeNew_z(k,i,:)=grayCode_z(k,b,:);  
    grayCodeNew_z(k,i+M/2,:)=grayCode_z(k,b,:);
    end
    y1(1,b)=0; %assign 0 to columns specified by 'b'in first row of y1
end

if same_num > 15   %fitness function
    for i=3:1:N
    x1(:,i)=rand(1,AP_num)*m;    
    x2(:,i)=uint8(x1(:,i)/m*PT);
    z1(:,i) = rand(1,AP_num)*n;
    z2(:,i)=uint8(z1(:,i)/n*PT);
    for j = 1:AP_num
    grayCodeNew_x(j,i,:)=num2gray(x2(j,i),L);
    grayCodeNew_z(j,i,:)=num2gray(z2(j,i),L);
    end  
    end
    same_num = 0;
end
    for i=2:1:M/2  %Crossover 
        p=unidrnd(L);  %uniform discrete random variables with maximum 8
        if rand()<Pc        
           
            for j=p:1:L   %modify bits if rand () reaches crossover probability
                for k =1:AP_num   %for each new AP
                temp=grayCodeNew_x(k,i,j);     
                grayCodeNew_x(k,i,j)=grayCodeNew_x(k,M-i+1,j); 
                grayCodeNew_x(k,M-i+1,j)=temp;  %new gray coded string stored in temp
                temp_z=grayCodeNew_z(k,i,j);   %similarly for z
                grayCodeNew_z(k,i,j)=grayCodeNew_z(k,M-i+1,j);
                grayCodeNew_z(k,M-i+1,j)=temp_z;
                end
            end
        end
    end
    %%
    for i=2:1:M  
        for j=1:1:L
            if rand()<Pm  %Mutation
                for r = 1:AP_num   
                grayCodeNew_x(r,i,j)=dec2bin(1-bin2dec(grayCodeNew_x(r,i,j)));    
                grayCodeNew_z(r,i,j)=dec2bin(1-bin2dec(grayCodeNew_z(r,i,j)));
                end
            end
        end
    end
%%
    for i=1:1:M
        for r = 1:AP_num
        x4(r,i)=gray2num(grayCodeNew_x(r,i,:));   %bin to dec conversion for gdop calc in Y3
        z4(r,i)=gray2num(grayCodeNew_z(r,i,:));
        end
    end
    x3=double(x4).*m/PT;    %divide by total population
    z3=double(z4).*n/PT;        
    y3 = value_calc_GDOP_0601(x3,z3,N,AP_num,RP_loc,AP_init_num);
 
[E,F] = min(y3);  %new AP's coordinate specified by F and determined through E which is minimum value of y3
plot(t,E,'b*')  %first figure plots GDOP 'E' against t 
hold on  
G_E(t+1)= E;  %updates for increments in t 
if G_E(t+1) == G_E(t)  
    same_num = same_num +1;  
else
    same_num = 0;  
end

    for i=1:1:N   %for 60 RPs   
        [a,b]=min(y3);  
        x1(:,i)=x3(:,b);  %new AP coordinates with reference to RPs' coordinates
        grayCode_x(:,i,:)=grayCodeNew_x(:,b,:);  
        z1(:,i)=z3(:,b);
        grayCode_z(:,i,:)=grayCodeNew_z(:,b,:);
        y3(1,b)=0;  
    end
end
 x1;  
 z1;
y1 =  value_calc_GDOP_0601(x1,z1,N,AP_num,RP_loc,AP_init_num);  %GDOP at new AP location
[E,F] = min(y1);   %minumum of y1 is the GDOP of new AP coordinates specified by F
x1(:,F) %x-coordinate for new AP 
z1(:,F) %z-coordinate for new AP
E       %GDOP

figure(2)
plot(x1(:,F),z1(:,F),'r*'); %plot new AP
hold on
plot(2,5,'*');  %plot existing AP 1
hold on
plot(14,5,'*')   %existing AP 2
hold on
plot(8,25,'*')   %existing AP 3
hold on
 

%% localize with the deployed AP
test_num = 110*35;
m = 22;
n =7;
a(1,1)  = 0.5;
a(1,2)  = 0.5;
%%%%%%%
power_trans = 24;
nh = 2.3;
M = 2;
loss_a = 3.1;
N = 2;
loss_b = 3.1;
for i =1:m
    for j = 1:3
        k = (j-1)*22 +i;
        a(k,1) = 0.5 + 2*(i-1);
        a(k,2) = 0.5 + 2*(j-1);
    end
end

for i = 1:22
    for j = 1:2
        p = (2+j)*22 + i;
        q = 2*22 + i;
        a(p,1) = a(q,1);
        a(p,2) = a(q,2) + 1*j;
    end
end
for i = 1:22
    for j = 1:2
        p = 4*22 + i;
        q = (4+j)*22 + i;
        a(q,1) = a(p,1);
        a(q,2) = a(p,2) + 2*j;
    end
end
RP_loc = a; 
temp_a=[1:14];
temp_b=[1:25];
a=mean(temp_a);
b=mean(temp_b);
    new_AP_loc=[min(temp_a),b;max(temp_a),b];
pos_a = new_AP_loc(1,:);
pos_c = new_AP_loc(2,:);
    new_AP_loc_l=[a,max(temp_b);a,min(temp_b)];
pos_b = new_AP_loc_l(1,:);
pos_d = new_AP_loc_l(2,:);

AP_loc = [pos_a(:,1) pos_c(:,1) pos_b(:,1) pos_d(:,1);pos_a(:,2) pos_c(:,2) pos_b(:,2) pos_d(:,2)];  
for i = 1:AP_total
    for j =1: 154
        dis = sqrt((RP_loc(j,1) - AP_loc(1,i))^2 + (RP_loc(j,2) - AP_loc(2,i))^2);
        path_loss(i,j) = MWF_model(M,loss_a,N,loss_b,dis,nh);
        power_rec(i,j) = power_trans - path_loss(i,j) + lognrnd(0.1,1);
    end
end
for i = 1:test_num
    l_y = floor(i/110);
    l_x = mod(i,35);
    pos_x = 0.4*l_x;
    pos_y = 0.4*l_y;    
               
    local_real(1) = pos_x;
    local_real(2) = pos_y;
    for j = 1:AP_total
        dis_real = sqrt((pos_x - AP_loc(1,j))^2 + (pos_y - AP_loc(2,j))^2);
        path_loss_real(j) = MWF_model(M,loss_a,N,loss_b,dis_real,nh);
        power_rec_real(j) = power_trans - path_loss_real(j) + lognrnd(0.1,1);
    end
  
    K =3;
    local_error(i) = knn_algorithm_2212(RP_loc,power_rec,power_rec_real,local_real,K,AP_total);

end

local_error(i)

figure(3)
plot(AP_loc(1,1),AP_loc(2,1),'b*')
hold on
plot(AP_loc(1,2),AP_loc(2,2),'b*')
hold on
plot(AP_loc(1,3),AP_loc(2,3),'b*')
hold on
plot(AP_loc(1,4),AP_loc(2,4),'r*')
hold on
axis([1 14 1 25]);


H1 = figure(4);
cdfplot(local_error)
saveas(H1,'position_error');
