function [TT_updated min_cond_vec] = nsht_ordered_theta(L)
%need to alter so reads in the .mat file corresponding to each sampling
%scheme

% nsht_ordered_theta - Determine the optimal placement of rings along theta
%
% Either reads in existing .mat files with optimal locations or performs
% step 4 of optimal placement method to create a .mat file with optimal
% locations in it
%
% Default usage is given by
%
%   f = nsht_ordered_theta(L)
%
% where L is the harmonic band-limit. TT_updated denotes the \theta vector 
% contains sampling positions along \theta. min_cond_vector contains the condition number
% of the matrix Y_m. See the following paper for details.
%
%       
% Original Author: Zubair Khalid
% Altered by: Alice Bates, July 2014
% NSHT package to perform spherical harmonic transforms
%


% Check arguments.
% if ~isreal(L)
%       error('Harmonic band-limit must be real');
% end 


try
    
    %load(['theta_locations/theta_step4_MAX_ANGLE_pi_L_',num2str(L),'.mat']);
    load(['theta_locations/theta_antipodal_NEW_L_',num2str(L),'.mat']);
  
    return
catch err
     
    %% carry out step 4 of optimal placement method
    %load in theta locations created from carrying out first three steps of
    %optimal placement method
    load(['theta_locations/theta_ring_locations_step1_3_MAX_ANGLE_pi_L',num2str(L),'.mat']);
    TT = theta_final;  
    [~,temp_index ] = min(abs(theta_final - 0.5*pi)); %find location of the largest ring (closest to pi/2)
    TT_ordered_index = [temp_index];
    min_cond_vec = [];
    %%


    for m = L-3:-2:0
            
             [P_mat_1 Sc_1] = nsht_legmat_mex(TT, L, m); 
             [P_mat_2 Sc_2] = nsht_legmat_mex(TT, L, m-1); %m odd            
            index_1 = 1:2:(L-m);
            index_2 = 2:2:L-m;
            
            Y_mat_1 = zeros(length(index_1),length(index_1));
            Y_mat_2 = zeros(length(index_2),length(index_2));

            PP_1 = 10.^Sc_1.*P_mat_1;
            PP_1 = PP_1(index_1,:); %keep only even degree rows

            PP_2 = 10.^Sc_2.*P_mat_2;
            PP_2 = PP_2(index_2,:); %keep only even degree rows
            
            
            for jj=1:1:length(TT_ordered_index)
                Y_mat_1(:,jj) = PP_1(:,TT_ordered_index(jj));
            end

            for jj=1:1:length(TT_ordered_index)
                Y_mat_2(:,jj) = PP_2(:,TT_ordered_index(jj));
            end

            ii_vec=[];
            cond_vec_1 = [];
            cond_vec_2 = [];
            for ii=1:1:length(TT)


                if isempty(find(TT_ordered_index==ii))
                    Y_mat_1(:,end) = PP_1(:,ii);
                    Y_mat_2(:,end) = PP_2(:,ii);
                    ii_vec = [ii_vec ii];
                    cond_vec_1 = [cond_vec_1 cond(Y_mat_1)];
                    cond_vec_2 = [cond_vec_2 cond(Y_mat_2)];
                end
            end

            [M,II] =min(cond_vec_1 + cond_vec_2);

            min_cond_vec = [min_cond_vec M];
            TT_ordered_index = [ ii_vec(II) TT_ordered_index];

    end

    TT_updated = TT(TT_ordered_index);
    
    save(['theta_locations/theta_antipodal_L_' num2str(L) '.mat'], 'TT_updated', 'min_cond_vec');
    %%save theta locations - NOTE MAY HAVE TO CHANGE VALUE OF MAX ANGLE
   % save(['theta_locations/theta_step4_MAX_ANGLE_pi_L_' num2str(L) '.mat'], 'TT_updated', 'min_cond_vec');
   
end %% end try/catch

end


