
clear;
%2025.8.14，从plink格式转换成.vcf格式，样本名称重复，改名。

path='';
fname='MouseW5B_v7.vcf';
fname='MouseW5B_v7.vcf';
fname='mgp_REL2021_v8.vcf';
fname='mgp_REL2021_v8hSNP.vcf';
fname='mgp_REL2021_v8X46.vcf';
fname='mgp_REL2021_v8hSNPX46.vcf';
fname='mgp_REL2021_46_v7hSNP.vcf';
fname='mgp_REL2021_46_v7hSNPX46JKS.vcf';
fname='mgp_REL2021_46_LDp4v7hSNP.vcf';
 fname='mgp_REL2021_46_LDp4v7punch256.vcf';
 fname='mgp_REL2021_LdW1p4v7.vcf';

fwname=[fname(1:end-4),'R.vcf'];

% 打开文件
    t=clock;
    fprintf('start file %s at time %d:%d:%2.1f\n',fname,t(4),t(5),t(6));
    fn=[path,fname];
    fp = fopen(fn,'r');
    
    data_vcf=fread(fp,'char');
    %数据文件分行；
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
    Lidx=Lidx-1;
    %计算文件中每行的样本数量
    snpv=zeros(LineAll,Nsample);
    snpp=zeros(1,LineAll);
    chrom=snpp;
    t10=10.^[0:20];
    snp_idx=0;
    Nhsnp=0;

    fwp = fopen([path,fwname],'w');
    fwrite(fwp,data_vcf(1:Line(Lidx-1)));
    %plink转换过来的样本名称有重复，需要改名；
    line=data_vcf(Line(Lidx-1)+1:Line(Lidx));
    sample=find(line==9);
    lineA=line(1:sample(8));
    Ns=length(sample)-9;
    for n=1:Ns,
        name=line(sample(n+8):sample(n+9)-1);
        ta=length(name);
        name=name(1:ta/2);
        lineA=[lineA',name']';
    end;
    name=line(sample(n+9):end-1);
    ta=length(name);
    name=name(1:ta/2);
    lineA=[lineA',name',10]';
    fwrite(fwp,lineA);

    fwrite(fwp,data_vcf(Line(Lidx)+1:end));
    fclose(fwp);
    return;
    
    
Nh=length(hgen_id);
   for Lidx=Lidx:length(Line)-1,
        Ldata=data_vcf(Line(Lidx)+1:Line(Lidx+1)); 
        List=find(Ldata==9);%find tab charictor
        tt=Ldata(List(9:end)+1);
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
        %已有SNP位置，查找SNP是否属于候选基因
        for n=1:Nh,
            if(chrom(snp_idx)==hgen_line(1,n))%比较基因号
                if(snpp(snp_idx)>hgen_line(2,n) && snpp(snp_idx)<hgen_line(3,n))%比较起始，结束位置
                    fwrite(fwp,Ldata);
                    Nhsnp=Nhsnp+1;
                end;
            end;
        end;
    end;
    fclose(fp);
    fclose(fwp);
%     %根据实际读取到的数据，截取变量大小
%     snpv=snpv(1:snp_idx,:);
%     snpp=snpp(1:snp_idx);
%     chrom=chrom(1:snp_idx);
% %     fsname=sprintf('W3core%d.mat',fc);
%     fsname=[path,fname(1:end-4),'.mat'];
%     save(fsname,'snpv','chrom','snpp','SampleName');
    t=clock;
    fprintf('%d varint have save to .vcf file %d:%d:%2.1f\n',Nhsnp,t(4),t(5),t(6));
% % end;
