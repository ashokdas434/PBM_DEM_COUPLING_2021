function dN = discrete_PBM(t, N, f_c, beta_star, psi, vol_range, avg_vol, PBM_trigger, t_lim)

%% beta formulation
NN = N*N';
beta_st_NN = beta_star .* NN;
beta_0 = (f_c * sum(N))/ sum(beta_st_NN(:));
beta = beta_star * beta_0 * psi;

%% Discrete PBE calculations
Birth = zeros(length(vol_range),1); Death = zeros(length(vol_range),1);

for i=1:length(vol_range)
    if i~= 1
        for j=1:(i-1)
            Birth(i) = Birth(i) + 0.5* beta(i-j,j)*N(i-j)*N(j);
        end
    end
    
    if i~= length(vol_range)
        for j=1:length(vol_range)-1
            Death(i) = Death(i) + beta(i,j)*N(i)*N(j);
        end
    end
end

dN = Birth - Death;
   
fprintf('t=%1.2f| N_p=%1.4f| Vol=%1.4f|\n',t,sum(N),vol_range*N)
end
