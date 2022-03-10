function [value,isterminal,direction] = event_function(t,N, f_c, beta_star, psi, vol_range, avg_vol, PBM_trigger,t_lim)

avg_V = (vol_range*N)/ sum(N); % avg vol after PBM
value = min(PBM_trigger *avg_vol - avg_V,t_lim-1e-5-t);  % when value = 0, an event is triggered
isterminal = 1; % terminate after the first event
direction = 0;  % get all the zeros
end
