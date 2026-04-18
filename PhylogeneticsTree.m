if(1)%0：手动中断，继续运行，1：从头开始筛选。
clear;
%NJ法进化树生成。输入为SNP，生成进化树，比较和文件上是否相同，可以在生成中检查是否和参考树一致，如果不一致就停下。
%2025.11.5，，

t=clock;
fprintf('begin processing %d:%d:%2.1f\n',t(4),t(5),t(6));

load('mgp_REL2021_LdW1p4v7R_L315669');
fmname='GenList_hSNPrefM39.mat';
load(fmname);

[N,M]=size(snpv);
fprintf('%d varint %d samples load from file\n',N,M);
No=N;

snp_punch=snpv;
Nh=length(hgen_id);

new_line=1;
%**********************手动参数配置区域***************************************
% 从30000基础上往下继续筛选。不判断是否是hSNP 
% load('X:\WorkData\MGItree\MGI46_CS_phy_L2999_11072033');%从昨天的结果继续往下筛。
% load('X:\WorkData\MGIphylogenetic\MGI_CS_phy_L29994_11121025');%从昨天的结果继续往下筛。
% load('X:\WorkData\MGIphylogenetic\MGI_CS_phy_L3000_11121025');%从昨天的结果继续往下筛。
Nresv=300;%保留点数
%*******************************************************************
if(Nresv>=1000)%筛选到1000点的时候，需要降低打孔数量。
    Npunch=4;
else
    Npunch=4;
end;
Pref=1;%第一次运行的点是基准点，不打孔。
Phy_Sort=zeros(2,M-2);

% while N>30000%先保留3万个，进一步删除要改变策略。
while N>Nresv%从上次的3万个继续往下筛选，保留3K，300个。
    if(Pref==0)%第一次运行，不打孔，作为基准。
        %随机删除n个点位，为了有点讲究，先不删除hSNP。等到剩余30k，再不分，随机删除
        pp=rand(1,Npunch)*(N-1)+1;
        pp=round(pp);
        %核对一下是否在hSNP里面。
        %从30000的基础上往下继续筛选，不判断是否是hSNP了。
        if(Nresv>=20000)    
            pNh=0;
            for np=1:Npunch,
                pChrom=chrom(pp(np));
                pSnp=snpp(pp(np));
                for nh=1:Nh,                
                    if(pChrom==hgen_line(1,nh) && pSnp>hgen_line(2,nh) && pSnp<hgen_line(3,nh))
                       pNh=pNh+1;
                       continue;
                    end;
                end;
            end;
            if(pNh>0)  
    %             fprintf('hSNP find\n');
                continue;    %有snp在hSNP集合里面，跳过重来。
            end;
        end;    
    %     校验通过，执行删除和测试。
        snp_punch(pp,:)=[];
     end;%if(Pref==0)，不打孔    
    
    %1、生成Pairwise Distances矩阵，和软件生成的做比对。进化树基于这个矩阵搞
    PDsnp=zeros(M,M);
    for m=2:M,
        for mm=1:m-1
            ta=snp_punch(:,m)-snp_punch(:,mm);
            tb=sum(abs(ta));
            PDsnp(m,mm)=tb/N;
            PDsnp(mm,m)=tb/N;
        end;
    end;
    %进化树构建循环，计算M矩阵，并记录最小值
    %逐步比对进化树构建，如果和基线不符合，停止迭代，重新打孔。
    phy_sort=zeros(2,M-2);
    r=sum(PDsnp);
    Phy_Ok=1;%默认构建树和基线一致，不一致时修改。
    Ma=M;
    for k=1:M-2,
        MM=zeros(Ma,Ma);
        Mmin=1;
        for m=2:Ma,
            for mm=1:m-1
                Mi=PDsnp(m,mm)-(r(m)+r(mm))/(Ma-2);
                MM(m,mm)=Mi;
                MM(mm,m)=Mi;
                if(Mi<Mmin)
                    m1=mm;
                    m2=m;
                    Mmin=Mi;
                end;
            end;
        end;
        phy_sort(1,k)=m1;
        phy_sort(2,k)=m2;
%         对比基线，看看变没变。
        if(Pref==0)
            if(Phy_Sort(1,k)~=phy_sort(1,k) || Phy_Sort(2,k)~=phy_sort(2,k))
                Phy_Ok=0;%进化树变形指示。
                break;
            end;
        else
            Phy_Sort(1,k)=m1;
            Phy_Sort(2,k)=m2;
%             Min_Dist_ref(k)=min_dist;
        end;

        %合并叶节点，刷新距离矩阵
        DUsnp=zeros(1,Ma);
        for m=1:Ma,
            DUsnp(m)=(PDsnp(m1,m)+PDsnp(m2,m)-PDsnp(m1,m2))/2;
        end;
        %m1小，m2大。DUsnp放到PDsnp(m1),PDsnp(m2)删除。
        PDsnp(m1,:)=DUsnp;
        PDsnp(:,m1)=DUsnp';
        PDsnp(m2,:)=[];
        PDsnp(:,m2)=[];
        
    %     记录m1,m2作为判断进化树是否正确的判据。
    %     name1=SampleName{m1};
    %     name2=SampleName{m2};
    %     SampleName{m1}=strcat(name1,name2);
    %     SampleName(m2)=[];
    %     fprintf('%s+%s\n',name1,name2);
    %     fprintf('%d+%d\n',m1,m2);
        r=sum(PDsnp);
        Ma=Ma-1;
    end;
    if(Pref==0)%判断是否是在构建发育树的基准
        if(Phy_Ok==1)%执行完构建发育树，中间结果都对
    % %     效果好，删了。
            chrom(pp)=[];
            snpp(pp)=[];
            snpALT(pp)=[];
            snpREF(pp)=[];   
            snpv=snp_punch;
             N=length(snpALT);
            fprintf('try %d times, %d SNPs left\n',new_line,N);        
             new_line=1;

        else
            %删除之后效果不好，补回来。
            snp_punch=snpv;
    %         fprintf('%d/45  ',k);
            new_line=new_line+1;
    %         if(new_line>16)
    %             new_line=0;
    %             fprintf('\n');
    %         end;        
        end;
    else
        Pref=0;
    end;       
end;
end;%if(0)手动中断运行情况下，保存结果。
fm_name=sprintf('MGI_CS_phy_L%d_%02d%02d%02d%02d',N,t(2),t(3),t(4),t(5));
save(['X:\WorkData\MGIphylogenetic\',fm_name]);%保存运行结果，可以继续运行。
% load('X:\WorkData\MGItree\MGI46_CS_phy_L29999_11061634');%仅用.mat生成.vcf文件
%生成.vcf文件
%生成文件内容，文件头手动加上去。
fv_name=sprintf('MGI_CS_phy_L%d_%02d%02d%02d%02d.vcf',N,t(2),t(3),t(4),t(5));
fn=[fv_name];
[Nsnp,Nsmp]=size(snpv);
fp=fopen(fn,'w');
fprintf(fp,'#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT');
for k=1:Nsmp,
    fprintf(fp,'\t%s',SampleName{k});
end;
fprintf(fp,'\n');

for k=1:Nsnp
    if(chrom(k)=='X')
        fprintf(fp,'X\t');
    else
         fprintf(fp,'%d\t',chrom(k));
    end;
    fprintf(fp,'%d\tC%dR%d\t%s\t%s\t',snpp(k),chrom(k),snpp(k),snpREF(k),snpALT(k));
    fprintf(fp,'.\tPASS\tTYPE=snp\tGT');
    for m=1:Nsmp
        ta=snpv(k,m);
        if(ta==0.5)
            fprintf(fp,'\t.');
        else
            fprintf(fp,'\t%d',ta);
        end;
    end;
    fprintf(fp,'\n');
end;
fclose(fp);
