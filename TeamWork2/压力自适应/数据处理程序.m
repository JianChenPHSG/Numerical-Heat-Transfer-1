clear;
clc;
close all;

FileName = '3-松弛因子-relax(1)=0.4';
Convergence_accuracy = 1E-5;

fidin=fopen([FileName,'.DAT']);
Line_Process=1;
DataProcess = zeros();
while ~feof(fidin)                                      % 判断是否为文件末尾               
   tline=fgetl(fidin);                                 % 从文件读行   
   size_tline=size(tline);
   if(size_tline(2)>100)
       tlinedouble = str2num(tline); %#ok<ST2NM>
       size_tlinedouble = size(tlinedouble);
       if(size_tlinedouble(2)==7)
           DataProcess(Line_Process,1)=round(tlinedouble(1,5)*10);
           DataProcess(Line_Process,2)=tlinedouble(1,7);
           Line_Process=Line_Process+1;
       else
           continue;
       end
   else
       continue;
   end
   if (mod(Line_Process,2500)==0)
       Line_Process %#ok<NOPTS>
   end
end
 
DataSize = size(DataProcess);
Relaxation_factor_number = max(DataProcess(:,1));
Convergence_Time = zeros(1,Relaxation_factor_number);
Convergence_Min = zeros(1,Relaxation_factor_number);
Data_Size_Of_Different_Relaxation_factor_number=zeros(Relaxation_factor_number+1,1);
for i=1:1:DataSize(1)
        Data_Size_Of_Different_Relaxation_factor_number(DataProcess(i,1)+1)=Data_Size_Of_Different_Relaxation_factor_number(DataProcess(i,1)+1)+1;
end

figure(1);
for i=1:1:Relaxation_factor_number
    Convergence_Time_Temp=find(DataProcess(sum(Data_Size_Of_Different_Relaxation_factor_number(1:i))+1:sum(Data_Size_Of_Different_Relaxation_factor_number(1:i+1)),2)<Convergence_accuracy, 1, 'first');
    Convergence_Min(i) = min(DataProcess(sum(Data_Size_Of_Different_Relaxation_factor_number(1:i))+1:sum(Data_Size_Of_Different_Relaxation_factor_number(1:i+1)),2));
    if(size(Convergence_Time_Temp)>=1)
        Convergence_Time(i)=Convergence_Time_Temp;
    else
        Convergence_Time(i) = inf;
    end
    loglog(DataProcess(sum(Data_Size_Of_Different_Relaxation_factor_number(1:i))+1:sum(Data_Size_Of_Different_Relaxation_factor_number(1:i+1)),2),'linewidth',1);
    grid on;
    Relaxation_factor = i/10;
    title(['Relaxation factor=',num2str(Relaxation_factor),' 收敛精度随迭代次数变化关系'],'FontSize',16,'Fontname', '黑体');
    xlabel('迭代次数','FontSize',12,'Fontname', '黑体');
    ylabel('收敛精度','FontSize',12,'Fontname', '黑体');
    set(gcf,'position',[100,100,1000,700]); 
    if ~exist(FileName,'dir')
    mkdir(FileName) % 若不存在，在当前目录中产生一个子目录
    end
    saveas(gcf,[FileName,'\','Relaxation factor=',num2str(Relaxation_factor),'.jpg']);
end
figure(2);
    i=1;
    loglog(DataProcess(sum(Data_Size_Of_Different_Relaxation_factor_number(1:i))+1:sum(Data_Size_Of_Different_Relaxation_factor_number(1:i+1)),2),'linewidth',1.5,'color','y');
    hold on;
    i=5;
    loglog(DataProcess(sum(Data_Size_Of_Different_Relaxation_factor_number(1:i))+1:sum(Data_Size_Of_Different_Relaxation_factor_number(1:i+1)),2),'linewidth',1.5,'color','m');
    i=6;
    loglog(DataProcess(sum(Data_Size_Of_Different_Relaxation_factor_number(1:i))+1:sum(Data_Size_Of_Different_Relaxation_factor_number(1:i+1)),2),'linewidth',1.5,'color','g');
    i=10;
    loglog(DataProcess(sum(Data_Size_Of_Different_Relaxation_factor_number(1:i))+1:sum(Data_Size_Of_Different_Relaxation_factor_number(1:i+1)),2),'linewidth',1.5,'color','b');
    title('收敛精度随迭代次数变化关系','FontSize',32,'Fontname', '黑体');
    xlabel('迭代次数','FontSize',24,'Fontname', '黑体');
    ylabel('收敛精度','FontSize',24,'Fontname', '黑体');
    set(gcf,'position',[100,100,2000,1400]); 
    hleg1 = legend('Relaxation factor= 0.1','Relaxation factor= 0.5','Self Adapting','Relaxation factor= 1');
    legend('boxoff')
    set(hleg1,'FontSize',24,'Fontname', 'Times New Roman','Location','NorthEast');
    grid on;
    saveas(gcf,[FileName,'\',[FileName,'Self Adapting'],'.tif']);
    hold off;

fprintf(['收敛精度的最小值为','\r']);
fprintf('RelaxationFactor\t收敛精度的最小值\r');
for i=1:1:Relaxation_factor_number
    Relaxation_factor = i/10;
    fprintf([num2str(Relaxation_factor),'\t',num2str(Convergence_Min(i)),'\r'])
end

fprintf(['收敛到精度小于',num2str(Convergence_accuracy),'时迭代次数与收敛因子的关系为','\r']);
fprintf('RelaxationFactor\t迭代次数\r');
for i=1:1:Relaxation_factor_number
    Relaxation_factor = i/10;
    fprintf([num2str(Relaxation_factor),'\t',int2str(Convergence_Time(i)),'\r'])
end

