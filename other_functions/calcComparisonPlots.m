function [no_output] = calcComparisonPlots(GPR,GPR_d)
% calcComparisonPlots creates a combined plot of the two GPR methods
% inputs must be two GPR structs
no_output = [];


% Main Comparison Figure (two row subplot)
figure
subplot(2,1,1)
errorbar(GPR.X_Star,GPR.Y_Star_Hat,GPR.Y_Star_Var,'ro'), hold on
plot(GPR.X,GPR.Y,'bo','MarkerFaceColor','b'), grid on, hold off
xlabel('X Values'), ylabel('Y Values')
legend('Prediction \mu and \sigma^2','Raw Data')
title(GPR.title_str)

subplot(2,1,2)
errorbar(GPR_d.X_Star,GPR_d.Y_Star_Hat,GPR_d.Y_Star_Var,'ro'), hold on
plot(GPR_d.X,GPR_d.Y,'bo','MarkerFaceColor','b'), grid on, hold off
xlabel('X Values'), ylabel('Y Values')
legend('Prediction \mu and \sigma^2','Raw Data')
title(GPR_d.title_str)


% Just Mean Values plot (one plot, mean data of both methods)
figure
plot(GPR_d.X_Star,GPR_d.Y_Star_Hat,'--or','MarkerFaceColor','r'), hold on
plot(GPR.X_Star,GPR.Y_Star_Hat,'--ok','MarkerFaceColor','k'), hold off, grid on
xlabel('X Values'), ylabel('Y Values')
legend('Approx GPR \mu','Full GPR \mu')
title('Full GPR vs. Approximate Method Mean Values')


end

