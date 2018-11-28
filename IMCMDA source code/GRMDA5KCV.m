clc,clear;
str = load(['.\1.miRNA-disease association\knowndiseasemirnainteraction.txt']);
[~,disease]=xlsread(['.\1.miRNA-disease association\disease name.xlsx']);
[~,miRNA]=xlsread(['.\1.miRNA-disease association\miRNA name.xlsx']);
% nd:the number of diseases
% nm:the number of miRNAs
% pp:the number of known diseae-miRNA associations
nd = max(str(:,2));
nm = max(str(:,1));
[pp,~] = size(str);
% FS:the functional similarity between m(i) and m(j)
% FSP:Functional similarity weighting matrix
% SS:the semantic similarity between d(i) and d(j).
% SSP:semantic similarity weighting matrix
%
FS = load(['.\4.miNA functional simialrity\functional similarity matrix.txt']);
FSP = load(['.\4.miNA functional simialrity\Functional similarity weighting matrix.txt']);             
SS1 = load(['.\2.disease semantic similarity 1\disease semantic similarity matrix 1.txt']);
SS2 = load(['.\3.disease semantic similarity 2\disease semantic similarity matrix 2.txt']);
SS = (SS1+SS2)/2;
SSP = load(['.\2.disease semantic similarity 1\disease semantic similarity weighting matrix1.txt']);
SSP_2 = zeros(383)+1;
index = find(0 == SS);
SSP_2(index) = 0;
%interaction: adajency matrix for the disease-miRNA association network
%interaction(i,j)=1 means miRNA j is related to disease i
interaction = zeros(nd,nm);
for i = 1:pp
    interaction(str(i,2),str(i,1)) = 1;
end
[kd_ori,km_ori] = gaussiansimilarity(interaction,nd,nm);                   %calculate gaussiansimilarity,返回处理完高斯核后的疾病（kd）和mirna(km)
[sd,sm] = integratedsimilarity(FS,FSP,SS,SSP_2,kd_ori,km_ori);               % Integrated similarity for diseases and miRNAs

interaction = interaction';


interaction_ori = interaction;
[score_ori] = GRMDA(interaction_ori,sd,sm);
F = score_ori;
index_1 = find(1 == interaction_ori);

 auc = zeros(1,100);
 for i = 1 : 100
    i
    indices = crossvalind('Kfold', pp, 5);%将数据样本随机分割为5部分
    for j = 1:5 %循环5次，分别取出第i部分作为测试样本，其余两部分作为训练样本
        
        index_2 = find(j == indices);
        interaction(index_1(index_2)) = 0;
    
        [kd,km] = gaussiansimilarity(interaction',nd,nm);
        [sd,sm] = integratedsimilarity(FS,FSP,SS,SSP_2,kd,km);
    
        [result]=GRMDA(interaction,sd,sm);
        score_ori(index_1(index_2)) = result(index_1(index_2));
        interaction = interaction_ori;
    end
    pre_label_score = score_ori(:);
    label_y = interaction_ori(:);
    auc(i) = roc_1(pre_label_score,label_y,'red');
 end
 auc_avg = mean(auc) 