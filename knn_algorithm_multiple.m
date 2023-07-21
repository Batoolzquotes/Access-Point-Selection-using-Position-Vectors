%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%KNN algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [local_error, local_ue] = knn_algorithm_multiple(RP_loc,power_rec,power_real,local_real,K,U,AP_num)
N = length(RP_loc(:,:,1)); 
temp = ones(U,2,N); 
local_total =zeros(U,1,2); 
for i = 1:N
    total_dis = 0;
    for j = 1:AP_num
    total_dis = total_dis + (power_rec(j,i) - power_real(j))^2;
    end
    for u = 1:U
    temp(u,2,i) = total_dis;  
    temp(u,1,1) = i;
    end
end
% end
for u = 1:U
temp(u)= temp(u)';
temp_a = sort(temp(u),2,'ascend');
for i = 1:K
    local_total(u,:,:) = local_total(u,:,:)+ RP_loc(u,temp_a(u,i,1),:);
end
end
local_ue = local_total/K;
local_error = sqrt((local_ue(1) - local_real(1))^2 + (local_ue(2) - local_real(2))^2);