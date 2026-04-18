import pandas as pd
from scipy.spatial.distance import pdist, squareform
import numpy as np

def read_custom_vcf(file_path):
    """
    读取自定义格式的VCF文件，彻底处理缺失值和格式问题
    :param file_path: 文件路径
    :return: 样本名称列表，基因型矩阵DataFrame（行为SNP，列为样本）
    """
    # 跳过开头的非数据行
    skip_rows = 0
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('##') or line.startswith('#CHROM'):
                skip_rows += 1
            else:
                break
    
    # 读取数据，前9列是基础信息（#CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT），第10列开始是样本基因型
    df = pd.read_csv(file_path, sep='\t', skiprows=skip_rows, header=None)
    
    # 修正：提取样本名称（第10列开始的列名，对应df的第9列索引及之后）
    # 注意：pd.read_csv(header=None)时，df.columns是数字索引。我们需要用df.iloc[0]来获取表头行，但这里已跳过#CHROM行。
    # 更稳健的做法：如果跳过了#CHROM行，则第一行数据就是第一个位点，样本名在df的第一行（索引0）的第9列之后。
    # 但根据标准VCF，跳过的#CHROM行就是表头，样本名在那行。既然我们用header=None读取，样本名就在第一行数据（df.iloc[0]）的第9列及之后。
    # 然而，我们的skiprows已经跳过了包含#CHROM的行，所以读取进来的df第一行就是第一个变异位点，没有样本名了。
    # 因此，我们需要在跳过行时，特别保留#CHROM行作为表头。

    # 重新设计读取逻辑：先找到#CHROM行
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    header_line_index = None
    for i, line in enumerate(lines):
        if line.startswith('#CHROM'):
            header_line_index = i
            break
    
    if header_line_index is None:
        raise ValueError("在VCF文件中未找到#CHROM行")
    
    # 使用#CHROM行作为表头读取数据
    df = pd.read_csv(file_path, sep='\t', skiprows=header_line_index, header=0)
    
    # 此时，df的列名就是：['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '样本1', '样本2', ...]
    # 提取样本名称（从第10列开始）
    sample_names = df.columns[9:].tolist()
    
    # 提取基因型数据（所有行，从第10列开始）
    genotype_data = df.iloc[:, 9:]
    
    # 处理缺失值：将"."替换为-1，同时处理可能存在的空格等无效字符
    genotype_data = genotype_data.replace('.', -1)
    genotype_data = genotype_data.replace(' ', -1)
    
    # 转换为数值类型，处理转换错误
    genotype_matrix = genotype_data.apply(pd.to_numeric, errors='coerce')
    
    # 填充NaN值为-1（缺失值标记），并转换为整数
    genotype_matrix = genotype_matrix.fillna(-1).astype(int)
    
    # 验证维度一致性：样本名称数量应等于基因型矩阵的列数
    if len(sample_names) != genotype_matrix.shape[1]:
        raise ValueError(f"维度不匹配！样本名称数量{len(sample_names)}，基因型矩阵列数{genotype_matrix.shape[1]}")
    
    print(f"成功读取 {len(sample_names)} 个样本，{genotype_matrix.shape[0]} 个SNP位点")
    print(f"缺失值数量: {(genotype_matrix.values == -1).sum()}")
    
    return sample_names, genotype_matrix


def calculate_hamming_distance(genotype_matrix):
    """
    计算汉明距离矩阵，正确处理缺失值
    :param genotype_matrix: 基因型矩阵DataFrame（行为SNP，列为样本）
    :return: 样本间的汉明距离矩阵
    """
    # 转换为numpy数组，并转置，使行为样本，列为SNP
    matrix = genotype_matrix.values.T  # 转置，现在 matrix.shape = (n_samples, n_snps)
    
    n_samples = matrix.shape[0]
    n_snps = matrix.shape[1]
    distance_matrix = np.zeros((n_samples, n_samples))
    
    print(f"开始计算汉明距离矩阵...")
    print(f"样本数: {n_samples}, SNP数: {n_snps}")
    
    for i in range(n_samples):
        for j in range(i, n_samples):
            # 获取两个样本的基因型
            seq1 = matrix[i]
            seq2 = matrix[j]
            
            # 创建有效位点掩码（两个位点都非缺失值）
            valid_mask = (seq1 != -1) & (seq2 != -1)
            
            # 计算有效位点数
            n_valid = np.sum(valid_mask)
            
            if n_valid > 0:
                # 提取有效位点
                valid_seq1 = seq1[valid_mask]
                valid_seq2 = seq2[valid_mask]
                
                # 计算汉明距离
                diff = np.sum(valid_seq1 != valid_seq2)
                distance_matrix[i, j] = diff
                distance_matrix[j, i] = diff
            else:
                # 如果没有有效位点，距离设为NaN
                distance_matrix[i, j] = np.nan
                distance_matrix[j, i] = np.nan
    
    return distance_matrix

def save_distance_matrix(distance_matrix, sample_names, output_path):
    """
    保存汉明距离矩阵为CSV文件
    :param distance_matrix: 汉明距离矩阵
    :param sample_names: 样本名称列表
    :param output_path: 输出文件路径
    """
    print("DEBUG: save_distance_matrix函数被调用")  # 添加调试信息
    
    # 验证维度一致性
    if distance_matrix.shape[0] != len(sample_names):
        raise ValueError(f"维度不匹配！距离矩阵行数{distance_matrix.shape[0]}，样本名称数量{len(sample_names)}")
    
    # 创建DataFrame
    df_distance = pd.DataFrame(distance_matrix, index=sample_names, columns=sample_names)
    
    # 保存到CSV
    df_distance.to_csv(output_path, index=True)
    
    print(f"汉明距离矩阵已保存到: {output_path}")
    print(f"矩阵形状: {df_distance.shape}")
    
    # 显示前5行5列的预览
    print("\n距离矩阵预览（前5行5列）:")
    preview = df_distance.iloc[:5, :5]
    print(preview)


def main():
    # 配置文件路径和输出路径
    input_file = "MGI_CS_phy_L297_03211715.vcf"
    output_file = "297/hamming_distance_matrix"
    input_file = "mgp_REL2021_LdW1p4v7R_L315669.csv"
    output_file = "315669/hamming_distance_matrix"
    input_file = "MGI_CS_phy_L512_11121025.vcf"
    output_file = "512/hamming_distance_matrix.csv"
    
    try:
        # 读取自定义VCF文件
        print("正在读取文件...")
        sample_names, genotype_matrix = read_custom_vcf(input_file)
        
        # 计算汉明距离矩阵
        print("正在计算汉明距离矩阵...")
        distance_matrix = calculate_hamming_distance(genotype_matrix)
        
        # 保存结果
        print("正在保存结果...")
        save_distance_matrix(distance_matrix, sample_names, output_file)
        
        print("处理完成！")
        
    except FileNotFoundError:
        print(f"错误：找不到文件 {input_file}")
    except Exception as e:
        print(f"处理过程中发生错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
