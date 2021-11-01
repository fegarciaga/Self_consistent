function [n_up, n_dn, n_updn]=measure_bp(Phi_T, phi_old, B_up, B_dn, Proj_k_half, t_bp, t_pop, N_up, N_par, N_sites)
  %% Implement backwards propagation with the stored values
  % calculate backpropagated bra
  Phi_left=Phi_T;
    for ii=1:t_bp+1
        Phi_left=Proj_k_half*Phi_left;
        for jj=1:N_sites
            Phi_left(jj,1:N_up)=Phi_left(jj,1:N_up)*B_up(t_bp-ii+2,jj);
            Phi_left(jj,N_up+1:N_par)=Phi_left(jj,1+N_up:N_par)*B_dn(t_bp-ii+2,jj);
        end
        Phi_left=Proj_k_half*Phi_left;
        if mod(ii,t_pop)==0
            [Phi_left(:,1:N_up), Trsh]=qr(Phi_left(:,1:N_up),0);
            [Phi_left(:,1+N_up:N_par), Trsh]=qr(Phi_left(:,1+N_up:N_par),0);
        end
    end
    %% Once the bra is back propagated, the observable is calculated the usual way
    % calculate inverse matrices in order to calculate Green's function
    %% Check when it doesn't work 
    inv_O_up=inv(Phi_left(:,1:N_up)'*phi_old(:,1:N_up));
    inv_O_dn=inv(Phi_left(:,1+N_up:N_par)'*phi_old(:,1+N_up:N_par));
    temp_up=phi_old(:,1:N_up)*inv_O_up;
    temp_dn=phi_old(:,N_up+1:N_par)*inv_O_dn;
    G_up=temp_up*Phi_left(:,1:N_up)';
    G_dn=temp_dn*Phi_left(:,N_up+1:N_par)';
    n_up=diag(G_up);
    n_dn=diag(G_dn);
    n_updn=(diag(G_up)).'*diag(G_dn);
end