
clear;
%2025.8.7，hSNP所在基因列表文件41588_2018_223_MOESM5_ESM.xlsx，
%，这个文件中的基因位置是基于参考序列M38，更新到M39的位置，用基因名称做索引。

%打开GENCODE M37中，基因名称和基因位置的映射关系文件
path='';
fname='gencode.vM37.annotation.mat';
load([path,fname]);
%载入hSNP所在基因列表文件
fnxls=[path,'41588_2018_223_MOESM5_ESM.xlsx'];
[data,text,raw]=xlsread(fnxls,'Sheet1');
[N,M]=size(text);
gL=length(gen_id);
fprintf('%d genes in xls file, %d gens in GENCODE file\n',N-5,gL);

hgen_line=zeros(3,N-5);
hgen_id=cell(1,N-5);
gen_idx=0;
for n=1:N-5,
    %读取xls文件基因名称，获取编号序列
    xls_gen_id=text{n+3,1};
    xls_gen_id=xls_gen_id(8:end);
    xls_gen_name=text{n+3,2};
    
    for l=1:gL,
        gen_ref=gen_id{l};
        gen_ref=gen_ref(1:11);
        if(strcmp(xls_gen_id,gen_ref))
            break;
        end;
        gen_ref=gen_name{l};
         if(strcmp(xls_gen_name,gen_ref))
            break;
        end;       
    end;
    if(l~=gL)
        %fprintf('%d trans success\n',n);
        gen_idx=gen_idx+1;
        hgen_line(1,gen_idx)=gen_line(1,l);
        hgen_line(2,gen_idx)=gen_line(2,l);
        hgen_line(3,gen_idx)=gen_line(3,l);     
        hgen_id{gen_idx}=gen_id{l};
    else
        fprintf('%d-%s trans failar\n',n,xls_gen_name);
    end;
    
end;
hgen_line=hgen_line(:,1:gen_idx);
hgen_id=hgen_id(1:gen_idx);
fprintf('%d gene %d trans success\n',N-5,gen_idx)
fsname='GenList_hSNPrefM39.mat';
save([path,fsname],'hgen_line','hgen_id');
