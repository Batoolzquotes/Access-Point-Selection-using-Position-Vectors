function A = GDOP_calc_wifi_AP_loc(RP_loc,AP_loc_calc,N,AP_num,nh) 
G_temp = 10*log10(2.71)*nh; %path loss to be multiplied 
for i = 1:N
   for j = 1:AP_num
       r(i,j) = sqrt((RP_loc(i,1) - AP_loc_calc(1,j))^2 + (RP_loc(i,2) - AP_loc_calc(2,j))^2);  %Euclidean distance between each RP and AP
       K(j,1) = -G_temp *(RP_loc(i,1) - AP_loc_calc(1,j))/(r(i,j))^2;   %H-matrix for DOP %%two columns four rows  60 matrices i.e DOP at each RP
       K(j,2) = -G_temp *(RP_loc(i,2) - AP_loc_calc(2,j))/(r(i,j))^2; 
   end
 %DOP equation  
B = K'*K;
B_temp = inv(B);
A(i) = sqrt(trace(B_temp));
end

