
clear;
%2025.8.7，读取基因标注文件gff3转成.mat文件加快后续处理速度。
Lread=1e7;
path='';
fname='gencode.vM37.annotation.gff3';

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

    while strcmp(line(1),'#')
        line=sprintf('%s ',data_vcf(Line(Lidx)+1:Line(Lidx+1)));
        Lidx=Lidx+1;
    end;
    Lidx=Lidx-1;

    %计算文件中每行的样本数量
    Fstart=ftell(fp);
    fseek(fp,0,'eof');
    Fend=ftell(fp);
    fseek(fp,Fstart,'bof');
    
    LineAll=Nline-Lidx+1;
    fprintf('this file has %d gens in it.\n',LineAll);
    gen_line=zeros(3,LineAll);
    gen_id=cell(1,LineAll);
    gen_name=cell(1,LineAll);
    t10=10.^(0:20);
    gen_idx=0;
    Sidx=0;


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
        if(Ldata(1)==35)
            continue;
        end;
        List=find(Ldata==9);%find tab charictor
        %only reserve gene, 
        chrom_type=Ldata(List(2)+1:List(3)-1);
        chrom_type=sprintf('%s',chrom_type);
        if(strcmp(chrom_type,'gene'))
            gen_idx=gen_idx+1;
            %chrom in 0~1 tab range 
            chrom_char=Ldata(4:List(1)-1)-48;
            if(length(chrom_char)==1)
                gen_line(1,gen_idx)=chrom_char;
            else
                gen_line(1,gen_idx)=chrom_char(1)*10+chrom_char(2);
            end;
            %chrom start / end position
            gstart=Ldata(List(3)+1:List(4)-1)-48;
            gend =Ldata(List(4)+1:List(5)-1)-48;
            Np=length(gstart);
            gen_line(2,gen_idx)=t10(Np:-1:1)*gstart;
            Np=length(gend);
            gen_line(3,gen_idx)=t10(Np:-1:1)*gend;
            %chrom name, ID
            gid=Ldata(List(8)+11:List(8)+23);
            gen_id{gen_idx}=sprintf('%s',gid); 
            %gene name
            ta=Ldata(List(8)+1:end);
            tb=find(ta==59);%按照分号分隔；
            tc=ta(tb(3)+1:tb(4)-1);%第三，第四分号之间是基因名称
             gen_name{gen_idx}=sprintf('%s',tc(11:end));%转字符串格式
            
        end;
    end;
    fclose(fp);

    %根据实际读取到的数据，截取变量大小
    gen_line=gen_line(:,1:gen_idx);
    gen_id=gen_id(1:gen_idx);

    fsname=fname(1:end-5);
    save([path,fsname,'.mat'],'gen_line','gen_id','gen_name');
    t=clock;
    fprintf('%d gene samples have save to .mat file %d:%d:%2.1f\n',gen_idx,t(4),t(5),t(6));

