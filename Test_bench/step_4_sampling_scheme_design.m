install_sht;

for L = 5:4:61
    TT = nsht_indexed_theta(40*L);

    TT_ordered_index = [40*L];
    min_cond_vec = [];

    for m = L-3:-2:0
            index_1 =1:2:(L-m);
            index_2 = 2:2:L-(m-1);
            
            [P_mat_1 Sc_1] = nsht_legmat_mex(TT, L, m); 
            
            Y_mat_1 = zeros(length(index_1),length(index_1));
            

            PP_1 = 10.^Sc_1.*P_mat_1;
            PP_1 = PP_1(index_1,:); %keep only even degree rows

            
            for jj=1:1:length(TT_ordered_index)
                Y_mat_1(:,jj) = PP_1(:,TT_ordered_index(jj));
            end

            
             
            if m~=0
                [P_mat_2 Sc_2] = nsht_legmat_mex(TT, L, m-1); %m odd            
                Y_mat_2 = zeros(length(index_2),length(index_2));
             
                PP_2 = 10.^Sc_2.*P_mat_2;
                PP_2 = PP_2(index_2,:); %keep only even degree rows
               for jj=1:1:length(TT_ordered_index)
                    Y_mat_2(:,jj) = PP_2(:,TT_ordered_index(jj));
                end
             end
            

            ii_vec=[];
            cond_vec_1 = [];
            cond_vec_2 = [];
            for ii=1:1:length(TT)

                if isempty(find(TT_ordered_index==ii))
                    Y_mat_1(:,end) = PP_1(:,ii);
                    ii_vec = [ii_vec ii];
                    cond_vec_1 = [cond_vec_1 cond(Y_mat_1)];
                     if m~=0
                         Y_mat_2(:,end) = PP_2(:,ii);
                        cond_vec_2 = [cond_vec_2 cond(Y_mat_2)];
                     end
                end
            end
            
            
            if m~=0
                [M,II] =min(cond_vec_1+ cond_vec_2);
            else
                [M,II] =min(cond_vec_1);
            end
            
            min_cond_vec = [min_cond_vec M];
            TT_ordered_index = [ ii_vec(II) TT_ordered_index];

    end

    TT_updated = TT(TT_ordered_index);
    
    save(['theta_locations\theta_antipodal_NEW_L_' num2str(L) '.mat'], 'TT_updated', 'min_cond_vec');
end