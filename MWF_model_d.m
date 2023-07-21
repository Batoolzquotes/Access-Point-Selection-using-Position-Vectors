%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Â·ËðÄ£ÐÍ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path_loss = MWF_model_d(M_l,loss_a,N_l,loss_b,d,nh)
L0 = 20;
path_loss = L0+M_l*loss_a + N_l*loss_b + 10*nh*log10(d);
end