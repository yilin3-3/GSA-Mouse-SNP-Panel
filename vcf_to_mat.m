
clear;
%2024.6，把.vcf文件转成.mat文件，保留SNP信息和位点位置信息。
%2025.10.12，增加读取REF和ALT类型，方便后续生成.vcf文件。

Lread=1e7;

path='';

fname='mgp_REL2021_v8.vcf';
fname='Jackson46NX.vcf';
fname='mgp_REL2021_v8hSNPX46R.vcf';
fname='mgp_REL2021_46_v7.vcf';
fname='mgp_REL2021_46_v7hSNPR.vcf';
fname='mgp_REL2021_46_LDp4v7R.vcf';
fname='MGI46_CSNoPunch.vcf';
 fname='MGI46_CS_phy_L257729_11111542.vcf';
 fname='mgp_REL2021_LdW1p4v7R_L315669.vcf';

% fname='mgp_REL2021_46_LDp4v7hSNPR.vcf';
% path='X:\WorkData\MouseW5B\';
% fname='MouseW5BR.vcf';
% fname='MouseW5BhSNP.vcf';

% fc=16;
% for fc=1:1
%     fn=sprintf('D:\\WorkData\\VcfW3All\\core%d.vcf',fc);
% 打开文件
    t=clock;
    fprintf('start file at time %d:%d:%2.1f\n',t(4),t(5),t(6));
    fn=[path,fname];
    fp = fopen(fn,'r');
%     MH_windows=300;%search MH with in 300bp
    lastread=0;
    data_vcf=fread(fp,Lread,'char');
    if(length(data_vcf)<Lread)
        lastread=1;
    end;  

    Lidx=1;
    Line=find(data_vcf==10);
    Nline=length(Line);

    line=sprintf('%s ',data_vcf(Line(Lidx)+1:Line(Lidx+1)));
    Lidx=Lidx+1;

    while strcmp(line(1:2),'##')
        line=sprintf('%s ',data_vcf(Line(Lidx)+1:Line(Lidx+1)));
        Lidx=Lidx+1;
    end;

    List=find(line==9);%find tab charictor
    %读取.vcf里面的样本名称行
    for l=1:length(List)-9;
        SampleName{l}=line(List(l+8)+1:List(l+9)-1);
    end;
    SampleName{l+1}=line(List(l+9)+1:end);
    Lidx=Lidx+1;
    LL=Line(Lidx+1)-Line(Lidx);
    Nsample=length(SampleName);
    Lfile=1e5;

    Fstart=ftell(fp);
    fseek(fp,0,'eof');
    Fend=ftell(fp);
    fseek(fp,Fstart,'bof');
    LineAll=round((Fend)/LL);
    fprintf('this file is about %d bytes,%d samples %d snp points\n',Fend,Nsample,LineAll);

    %计算文件中每行的样本数量

    snpv=zeros(LineAll,Nsample);
    snpp=zeros(1,LineAll);
    snpALT=zeros(1,LineAll);
    snpREF=zeros(1,LineAll);
    chrom=snpp;
    t10=10.^[0:20];
    snp_idx=0;
    Sidx=0;
Lidx=Lidx-1;

   while(1)
   
        if(Lidx==Nline)
            if(lastread==1)
                break;
            end;
            Fstart=ftell(fp);
            t=clock;
            fprintf('%2.2f%% have read %d:%d:%2.1f\n',Fstart*100/Fend,t(4),t(5),t(6));            
            ta=fread(fp,Lread,'char');
            data_vcf=[data_vcf(Line(Lidx)-1:end)',ta']';
            Line=find(data_vcf==10);
            Nline=length(Line);
            Lidx=1;
            if(length(ta)<Lread)
                lastread=1;
            end;

        end;
        Ldata=data_vcf(Line(Lidx)+1:Line(Lidx+1));
        Lidx=Lidx+1;    
        List=find(Ldata==9);%find tab charictor
        tt=Ldata(List(9:end)+1);
        %从字符转换成0，0.5，1；
        tt=tt-48;
        ta=find(tt==-2);
        tt(ta)=0.5;
        
        snp_idx=snp_idx+1;
        snpv(snp_idx,:)=tt;
        %转换SNP位置字段
        ta=Ldata(List(2)-1:-1:List(1)+1)-48;
        Np=length(ta);
        tb=ta'.*t10(1:Np);
        snpp(snp_idx)=sum(tb);
        %转换染色体编号
        ta=Ldata(List(1)-1:-1:1)-48;
        Np=length(ta);
        tb=ta'.*t10(1:Np);
        chrom(snp_idx)=sum(tb);        
        %读取ALT，REF
        snpREF(snp_idx)=Ldata(List(3)+1:List(4)-1);
        snpALT(snp_idx)=Ldata(List(4)+1:List(5)-1);
        
    end;
    fclose(fp);

    %根据实际读取到的数据，截取变量大小
    snpv=snpv(1:snp_idx,:);
    snpp=snpp(1:snp_idx);
    chrom=chrom(1:snp_idx);
    snpREF=snpREF(1:snp_idx);
    snpALT=snpALT(1:snp_idx);
%     fsname=sprintf('W3core%d.mat',fc);
    fsname=[path,fname(1:end-4),'.mat'];
    save(fsname,'snpv','chrom','snpp','SampleName','snpALT','snpREF');
    t=clock;
    fprintf('%d varint %d samples have save to .mat file %d:%d:%2.1f\n',snp_idx,Nsample,t(4),t(5),t(6));
% end;
PlotAllDistance(snpv,1);