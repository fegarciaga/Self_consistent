function [Phi_T, E_T,delta]=singleit(H_k,Den_dn,Den_up,U_eff,N_up,N_dn,N_par,N_sites)
    % Arrange the independent particle hamiltonian
    H_ip_up=H_k;
    H_ip_dn=H_k;
    % Note an effective interaction is used 
    for j=1:N_sites
        H_ip_up(j,j)=H_ip_up(j,j)+U_eff*Den_dn(j);
        H_ip_dn(j,j)=H_ip_dn(j,j)+U_eff*Den_up(j);
    end
    [psi_ip_up,E_ip_up] = eig(H_ip_up);
    E_ipv_up=diag(E_ip_up);
    [psi_ip_dn,E_ip_dn] = eig(H_ip_dn);
    E_ipv_dn=diag(E_ip_dn);
    Phi_T=horzcat(psi_ip_up(:,1:N_up),psi_ip_dn(:,1:N_dn));
    E_T=sum(E_ipv_up(1:N_up))+sum(E_ipv_dn(1:N_dn));
    n_r_up=diag(Phi_T(:,1:N_up)*(Phi_T(:,1:N_up))');
    n_r_dn=diag(Phi_T(:,N_up+1:N_par)*(Phi_T(:,N_up+1:N_par))');
    n_r_up=reshape(n_r_up,[1,N_sites]);
    n_r_dn=reshape(n_r_dn,[1,N_sites]);
    % Find U_eff which minimizes the following quantity
    delta=1/N_par*sqrt((n_r_up-Den_up)*(n_r_up-Den_up)'+(n_r_dn-Den_dn)*(n_r_dn-Den_dn)');
end