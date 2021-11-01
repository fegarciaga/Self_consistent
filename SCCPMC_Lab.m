function [E_ave,E_err,Den_up_ave,Den_dn_ave,savedFileName]=SCCPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,tx2,ty2,tz2,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,t_bp,t_pop,N_iter,N_max,suffix)
    N_sites=Lx*Ly*Lz;
    Den_up_ave=zeros(N_sites,1);
    Den_dn_ave=zeros(N_sites,1);
    N_updn=0;
    for i=1:N_max
        if i==1
            flag=1;
            flag_save=0;
        else
            flag=0;
            if i==N_max
                flag_save=1;
            else
                flag_save=0;
            end
        end
        [E_ave,E_err,Den_up_ave,Den_dn_ave,N_updn,savedFileName]=CPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,tx2,ty2,tz2,Den_up_ave,Den_dn_ave,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,t_bp,t_pop,flag,flag_save,N_iter,suffix);
    end
end