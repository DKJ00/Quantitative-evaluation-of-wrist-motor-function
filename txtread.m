%%
%Data = importdata("E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0105\04\0221.txt");

%namelist = dir ('E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0105\01\*.txt'    );

% path = 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据0105\01\';
% namelist = dir ([path,'*.txt']);
% l = length(namelist);
% %P = cell(1,l);%定义一个细胞数组，用于存放所有txt文件
% %namelist(i).name;%这里获得的只是该路径下的文件名，如1.txt是相对路径
% for i = 1:l
%     filename{i} = [path,namelist(i).name];%通过字符串拼接获得的就是绝对路径了
%     X{i}= load(filename{i});
% end
% Y = cell2mat(X');
% 
% %T = [namelist.name;Y] ;
% xlswrite('C:\Users\16939\Desktop\falcon Data\01.xls',Y);%(在matlab工作窗口中的数组);
% for o = 0:2        
%     for p = 0:2           
%         for q = 1:1
%                 i = i+1;
%                 g = [path,'0',num2str(o), num2str(p),num2str(q) '.txt']; 
%                 Z(i) = (load(g));%读入 .txt 文件置于 data 细胞中  
%         end
%     end
% end
%%

%data = cell(1,45); %建立细胞存储空间
Z = [];
path = 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据1209\02\';
namelist = dir ([path,'*.txt']);
l = length(namelist);

    path = 'E:\OneDrive - mail.ustc.edu.cn\LABVIEW\实验数据1209\02\';
   i=0;
for o = 0:2        
    for p = 0:2           
        
                i = i+1;
                g = [path,'0',num2str(o), num2str(p) '.txt']; 
                Z(i) = (load(g));%读入 .txt 文件置于 data 细胞中  
       
    end
end
Y = Z';
V = [];
V = Y(1:2:end);

%xlswrite('C:\Users\16939\Desktop\falcon Data\03.xls',Y);


