
import pandas as pd
import numpy as np

# 读取CSV文件
#dataFile = '315669/hamming_distance_ordered.csv'
# dataFile = '512/hamming_distance_ordered.csv'
dataFile = '297/hamming_distance_ordered.csv'
dataCSV = pd.read_csv(dataFile)

# 将数据转换为numpy数组
matrix = dataCSV.values

# 将上三角部分（不包括对角线）置为0
matrix = np.tril(matrix)  # 保留下三角部分（包括对角线）

header = ["", "C58_J", "C57L_J", "C57BR_cdJ", "C57BL_6NJ", "C57BL_10J", "C57BL_10SnJ", 
          "B10.RIII", "LEWES_EiJ", "ZALENDE_EiJ", "WSB_EiJ", "CAST_EiJ", "MOLF_EiJ", 
          "JF1_MsJ", "PWK_PhJ", "CZECHII_EiJ", "BTBR_T+_Itpr3tf_J", "LP_J", "129S5SvEvBrd", 
          "129S1_SvImJ", "129P2_OlaHsd", "KK_HiJ", "NZW_LacJ", "NZO_HlLtJ", "NZB_B1NJ", 
          "I_LnJ", "SM_J", "LG_J", "RF_J", "AKR_J", "ST_bJ", "BUB_BnJ", "NOD_ShiLtJ", 
          "FVB_NJ", "SWR_J", "SJL_J", "CE_J", "RIIIS_J", "MAMy_J", "PL_J", "NON_LtJ", 
          "QSi5", "QSi3", "A_J", "SEA_GnJ", "BALB_cJ", "BALB_cByJ", "DBA_2J", "DBA_1J", 
          "CBA_J", "C3H_HeJ", "C3H_HeH"]
# 添加表头
df_with_header = pd.DataFrame(matrix, columns=header)

# 保存处理后的矩阵到新的CSV文件
#outputFile = '315669/hamming_distance.csv'
# outputFile = '512/hamming_distance.csv'
outputFile = '297/hamming_distance.csv'
df_with_header.to_csv(outputFile, index=False)


print("上三角部分已置零并保存至", outputFile)
