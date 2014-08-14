%testBench.m creates numerical accuracy plots
%author: Alice Bates
%date: July 2014

% Install SHT 
install_sht



MAX_Errors = zeros(1,14);
MEAN_Errors = zeros(1,14);

% 
% global conditionNum ;
% conditionNum = zeros(1,L);

%check L odd
% if mod(L,2)==0
%     L = L -1;
% end

index = 1; %5:4:61
for L = 9%9:4:61
    %for j =1:10 %repeat experiment 10 times for each L


        %% Generate random antipodal signal
        flmt = zeros(1,L^2);
        for l = 1:2:L
             flmt(l^2-2*l+2:l^2) =  rand(1,(2*l-1)) + 1i*rand(1,(2*l-1));
        end

        %% inverse transform - works  
        ft  = nsht_inverse(flmt,L);


        %% forward transform
        [flmr T_mat_inv] = nsht_forward(ft,L);


        %% Absolute error
        Error_max = max(abs(flmt-flmr));

        %%MEAN error
        Error_mean = sum(abs(flmt-flmr))/L^2;  
        
        %store average RMS and MAX error for plotting
        MAX_Errors(index) = MAX_Errors(index) +  Error_max*0.1;
        MEAN_Errors(index) = MEAN_Errors(index) + Error_mean*0.1;

  % end
    
    %index = index + 1;
end



%% plot max and mean error
x = 5:4:61;
figure();
semilogy(MAX_Errors,'.-','MarkerSize',14,'Color','k','linewidth',2);
hold on;
semilogy(MEAN_Errors,'-x','MarkerSize',10,'Color','k','linewidth',2);

set(gca,'FontSize',12) 
set(gca,'Xtick',1:15,'XTickLabel',{'8','','16','','24','','32','','40','','48','','56','','64'});
xlabel('Band-limit, L');
ylabel('E_{max} or E_{mean}');
legend('E_{max}','E_{mean}','Location','northwest');
