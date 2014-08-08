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
    
    load(['theta_locations/theta_step4_MAX_ANGLE_pi_L_',num2str(L),'.mat']);
    
  
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
            
             [P Sc] = nsht_legmat_mex(TT, L, m); 
                        
            index = 1:2:(L-m);
            Y_mat = zeros(length(index),length(index));

            PP = 10.^Sc.*P;
            PP = PP(index,:); %keep only even degree rows


            for jj=1:1:length(TT_ordered_index)
                Y_mat(:,jj) = PP(:,TT_ordered_index(jj));
            end


            ii_vec=[];
            cond_vec = [];
            
            for ii=1:1:length(TT)


                if isempty(find(TT_ordered_index==ii))
                    Y_mat(:,end) = PP(:,ii);
                    ii_vec = [ii_vec ii];
                    cond_vec = [cond_vec cond(Y_mat)];
                end
            end

            [M,II] =min(cond_vec);

            min_cond_vec = [min_cond_vec M];
            TT_ordered_index = [ ii_vec(II) TT_ordered_index];

    end

    TT_updated = TT(TT_ordered_index);
    %%save theta locations - NOTE MAY HAVE TO CHANGE VALUE OF MAX ANGLE
    save(['theta_locations/theta_step4_MAX_ANGLE_pi_L_' num2str(L) '.mat'], 'TT_updated', 'min_cond_vec');
end %% end try/catch

end


