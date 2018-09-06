% f1<-mnorm<-function(x){
%   y<-0.3135405*dnorm(x,8.359357,10.93816)+0.3278223*dnorm(x, 19.685839,36.93783)
%   +0.3586372*dnorm(x,31.938943 ,101.40636)}
clear clc;
close all

figure
hold on
x=0:0.1:100;
plot(x,f(x), '-', 'LineWidth', 2)
plot(x,normpdf(x,8.32,sqrt(10.8)))
plot(x,normpdf(x,19.65 ,sqrt(37.9)))
plot(x,normpdf(x,32.19 ,sqrt(99.8)))

%% defination
function y=f(x)
    y = 0.31*normpdf(x,8.32,sqrt(10.8))+0.34*normpdf(x,19.65 ,sqrt(37.9))+0.35*normpdf(x,32.19 ,sqrt(99.8));
end
