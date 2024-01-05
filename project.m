sensex_data = readtable('CSVForDate (2).csv')

fulldata = sensex_data(:,2);
fulldata = table2array(fulldata);
fulldata = fulldata';
size = length(fulldata);
disp(size);

opendata = fulldata(1240:1980);
n1=length(opendata)


closedata=sensex_data(:,2);
closedata=table2array(closedata);
closedata=closedata';

open_mean=mean(opendata);
disp(open_mean);
open_var=var(opendata);
disp(open_var);

figure
subplot(2,1,1);
hist(opendata);
title("Histogram plot of the series");
subplot(2,1,2);
plot(opendata);
title("Plot of the series");

figure
subplot(2,1,1);
autocorr(opendata);
title("ACF - Open Data");
subplot(2,1,2);
parcorr(opendata);
title("PACF - Open Data");

open1(1)= opendata(1);
 for i=2:n1
     open1(i)=opendata(i) - opendata(i-1);
end

figure
subplot(2,2,1)
plot(open1);
title("Plot of Open data - 1 Difference");
subplot(2,2,2)
plot(opendata);
title("Plot of Open data - 0 Difference");
subplot(2,2,3)
histogram(open1);
title("Histogram - 1 Difference");
subplot(2,2,4)
histogram(opendata);
title("Histogram - 0 Difference");

figure
subplot(2,2,1)
autocorr(open1);
title("ACF-  1 Difference");
subplot(2,2,2)
autocorr(opendata);
title("ACF - 0 Difference");
subplot(2,2,3)
parcorr(open1);
title("PACF- 1 Difference");
subplot(2,2,4)
parcorr(opendata);
title("PACF - 0 Difference");

ARIMA_Theoretical= arima(0,1,0);
[ARIMA_Theoretical1,~,LogLikelihood]= estimate(ARIMA_Theoretical, opendata');

resid= infer(ARIMA_Theoretical1, opendata');                                                
predict= opendata +resid';            
figure
subplot(1,2,1);
plot(opendata);
title("Original");
grid on
subplot(1,2,2);
plot(predict);
title("ARIMA (0,1,0)");
grid on

figure
histogram(resid);

[aic,bic]= aicbic(LogLikelihood,2,495);
Model_table = "arima(0,1,0)";
aic_table = aic;
bic_table = bic;

m=2;
for i = 1:2
    for j = 1:2
        for k=1:2
            arima_model(m-1) = arima(i,j,k);
            [~,~,LoglikehoodE] = estimate(arima_model(m-1),opendata','display','off');
            [aic_inter,bic_inter] = aicbic(LoglikehoodE,2,250);
            Model_table(m) = ['arima(',num2str(i),',',num2str(j),',',num2str(k),')'];
            aic_table(m) = aic_inter;
            bic_table(m) = bic_inter;
            m=m+1;
        end
    end
end


Comparision = table(Model_table',aic_table',bic_table');
disp(Comparision);

aicmin=aic_table(1);
aicindex = 1;
for i=1:length(aic_table)
    if(aic_table(i)<aicmin)
        aicindex=i;
        aicmin=aic_table(i);
    end
end

disp(aic_table);

disp("Lowest AIC value at:");

disp(aicindex);

BICMin=bic_table(1);
BICindex = 1;
for i=1:length(bic_table)
    if(bic_table(i)<BICMin)
        BICindex=i;
        BICMin=bic_table(i);
    end
end

disp(bic_table);

disp("Lowest BIC value at:");

disp(BICindex);

disp(Model_table);

disp(Model_table(BICindex));

arima_practical= arima(2,1,2);                               
[arima_prac1,~,LogLikelihood2]= estimate(arima_practical, opendata');

resid2= infer(arima_prac1, opendata');        
prediction2= opendata +resid2';

figure
subplot(1,2,1);
plot(opendata);
title("Original");
grid on
subplot(1,2,2);
plot(prediction2);
title("ARIMA (2,1,2)");
grid on

a = opendata'
b=full_opendata(1981:end,1)
Md=estimate(arima_prac1,a);
f=forecast(Md,577,'Y0',a);

figure
plot(b)
hold on
plot(f,'k--','LineWidth',1.5)
xlim([0,577])
title('Prediction Error')
legend('Observed','Forecast','Location','northwest')
hold off

pmse = mean((b-f).^2)

