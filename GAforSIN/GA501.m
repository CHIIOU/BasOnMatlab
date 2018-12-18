
clc;
clear all;
close all;
global BitLength
global boundsbegin
global boundsend
bounds=[0.6 25];
bounds2=[0.6 6.4];
precision=0.0001; 
boundsbegin=bounds(:,1);
boundsend=bounds(:,2);
%计算如果满足求解精度至少需要多长的染色体
BitLength=ceil(log2((boundsend-boundsbegin)' ./ precision)); 
pcrossover=0.90; 
pmutation=0.09; 
popsize=25; 
Generationnmax=12; 
%初始种群建立
population=round(rand(popsize,BitLength));
%计算适应度,返回适应度Fitvalue和累积概率cumsump
[Fitvalue,cumsump]=fitnessfun(population);  
Generation=1;
while Generation<Generationnmax+1
   for j=1:2:popsize 
      %选择操作
      seln=selection(population,cumsump);
      %交叉操作
      scro=crossover(population,seln,pcrossover);
      scnew(j,:)=scro(1,:);
      scnew(j+1,:)=scro(2,:);
      %变异操作
      smnew(j,:)=mutation(scnew(j,:),pmutation);
      smnew(j+1,:)=mutation(scnew(j+1,:),pmutation);
   end
   population=smnew;  %产生了新的种群
   %计算新种群的适应度   
   [Fitvalue,cumsump]=fitnessfun(population);
   %记录当前代最好的适应度和平均适应度
   [fmax,nmax]=max(Fitvalue);
   fmean=mean(Fitvalue);
   ymax(Generation)=fmax;
   ymean(Generation)=fmean;
   %记录当前代的最佳染色体个体
   x=transform2to10(population(nmax,:));
  xx=boundsbegin+x*(boundsend-boundsbegin)/(power((boundsend),BitLength)-1);
   xmax(Generation)=xx;
   Generation=Generation+1
end
Generation=Generation-1;
Bestpopulation=xx
Besttargetfunvalue=targetfun(xx)
figure(1);
hand1=plot(1:Generation,ymax);
set(hand1,'linestyle','-','linewidth',1.8,'marker','*','markersize',6)
hold on;
hand2=plot(1:Generation,ymean);
set(hand2,'color','r','linestyle','-','linewidth',1.8,...
'marker','h','markersize',6)
xlabel('进化次数');ylabel('最大/平均适应度');xlim([1 Generationnmax]);
legend('最大适应度','平均适应度');
box off;hold off;
