function local_error = knn_algorithm_comp_c(RP_loc,power_rec_c,power_real_c,temp_user,K,U,select_APs)
N = length(RP_loc(:,1));
temp = ones(2,N);
local_total =zeros(1,2);
for i = 1:N
    total_dis_c = 0;
    for j = 1:select_APs
    total_dis_c = total_dis_c + (power_rec_c(j,i) - power_real_c(j))^2;
    end
    temp(2,i) = total_dis_c;  %%%%%%%%%%%%%
    temp(1,i) = i;
end
temp= temp';
temp_a = sortrows(temp,2);
for i = 1:K
    local_total = local_total + RP_loc((temp_a(i,1)),:);
end
local_ue = local_total/K;
for u = 1:U
local_error = (sqrt((local_ue(1) - temp_user(u,1))^2 + (local_ue(2) - temp_user(u,2))^2))/10;
end