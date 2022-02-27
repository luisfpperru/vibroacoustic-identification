function plot_measures(data_filename)
    close all
    data = load(data_filename); 
    frequencia = data(:,1);
    NPS_IS = data(:,2);
    NPS_IU = data(:,3); 
    NPS_INN = data(:,4);
    sigma_IU = data(:,5);
    sigma_INN = data(:,6);
    eta_IU = data(:,7);
    eta_INN = data(:,8);
    MAC_IU = data(:,9);
    MAC_INN = data(:,10);
    figure, plot(frequencia,[NPS_IS,NPS_IU,NPS_INN]), 
    title('Comparison between Sound Power Level'), 
    legend('SPL IS','SPL IU','SPL INN',4)
    axis([0,max(frequencia),0,max([NPS_IS;NPS_IU;NPS_INN])])
    figure, plot(frequencia,[sigma_IU,sigma_INN]), title('Comparison of Pearson Correlation'),
    legend('IU','INN')
    figure, plot(frequencia,[eta_IU,eta_INN]), title('Comparison of Spearman Correlation'),
    legend('IU','INN')
    figure, plot(frequencia,[MAC_IU,MAC_INN]), title('Comparison of MAC Correlation'),
    legend('IU','INN')
end
