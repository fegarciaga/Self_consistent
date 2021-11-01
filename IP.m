function [Phi_T, E_T]=IP(H_k,flag,N_up,N_dn,Den_up,Den_dn,U,N_iter,N_sites)
    N_par=N_up+N_dn;
    if flag==1
        [psi_nonint,E_nonint_m] = eig(H_k);
        E_nonint_v=diag(E_nonint_m);
        % assemble the non-interacting single-particle orbitals into a Slater determinant:
        Phi_T=horzcat(psi_nonint(:,1:N_up),psi_nonint(:,1:N_dn));
        % the kinetic energy of the trial wave function
        E_K=sum(E_nonint_v(1:N_up))+sum(E_nonint_v(1:N_dn));
        % the potential energy of the trial wave function
        n_r_up=diag(Phi_T(:,1:N_up)*(Phi_T(:,1:N_up))');
        n_r_dn=diag(Phi_T(:,N_up+1:N_par)*(Phi_T(:,N_up+1:N_par))');
        E_V=U*n_r_up'*n_r_dn;
        % the total energy of the trial wave function = the initial trial energy
        E_T = E_K+E_V;
    else
        for i=1:N_iter
            if i==1
                U_eff_old2=U;
                U_eff_old1=U-0.05*U;
                [Phi_T, E_T,delta_old]=singleit(H_k,Den_dn,Den_up,U_eff_old2,N_up,N_dn,N_par,N_sites);
                [Phi_T, E_T,delta_new]=singleit(H_k,Den_dn,Den_up,U_eff_old1,N_up,N_dn,N_par,N_sites);
                ddelta=(delta_new-delta_old)/(U_eff_old1-U_eff_old2);
                U_eff_old2=U_eff_old1;
                U_eff_old1=U_eff_old1-ddelta;
                delta_old=delta_new;
            else
                [Phi_T, E_T,delta_new]=singleit(H_k,Den_dn,Den_up,U_eff_old1,N_up,N_dn,N_par,N_sites);
                ddelta=(delta_new-delta_old)/(U_eff_old1-U_eff_old2);
                U_eff_old2=U_eff_old1;
                U_eff_old1=U_eff_old1-ddelta;
                delta_old=delta_new;
            end
        end
    end
end