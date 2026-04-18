import pandas as pd

def reorder_hamming_matrix(csv_file, new_order):
    # 读取CSV文件，保留表头
    df = pd.read_csv(csv_file, index_col=0)
    
    # 检查新顺序中的品系是否都在原矩阵中
    missing_strains = [strain for strain in new_order if strain not in df.index.tolist() + df.columns.tolist()]
    if missing_strains:
        raise ValueError(f"以下品系不在原矩阵中：{missing_strains}")
    
    # 重新排序列（先处理列，再处理行）
    df = df[new_order]
    
    # 重新排序行
    df = df.loc[new_order]
    
    return df

if __name__ == "__main__":
    # 替换为你的CSV文件路径
    csv_file = "315669/hamming_distance_matrix.csv"
    csv_file = "512/hamming_distance_matrix.csv"
    csv_file = "297/hamming_distance_matrix.csv"
    output_csv='315669/hamming_distance_ordered.csv'
    output_csv='512/hamming_distance_ordered.csv'
    output_csv='297/hamming_distance_ordered.csv'

    # 目标品系顺序
    new_order = [
        "C58_J", "C57L_J", "C57BR_cdJ", "C57BL_6NJ", "C57BL_10J", "C57BL_10SnJ", "B10.RIII", "LEWES_EiJ", "ZALENDE_EiJ",
        "WSB_EiJ", "CAST_EiJ", "MOLF_EiJ", "JF1_MsJ", "PWK_PhJ", "CZECHII_EiJ", "BTBR_T+_Itpr3tf_J", "LP_J", "129S5SvEvBrd",
        "129S1_SvImJ", "129P2_OlaHsd", "KK_HiJ", "NZW_LacJ", "NZO_HlLtJ", "NZB_B1NJ", "I_LnJ", "SM_J", "LG_J", "RF_J",
        "AKR_J", "ST_bJ", "BUB_BnJ", "NOD_ShiLtJ", "FVB_NJ", "SWR_J", "SJL_J", "CE_J", "RIIIS_J", "MAMy_J", "PL_J",
        "NON_LtJ", "QSi5", "QSi3", "A_J", "SEA_GnJ", "BALB_cJ", "BALB_cByJ", "DBA_2J", "DBA_1J", "CBA_J", "C3H_HeJ",
        "C3H_HeH"
    ]
    
    # 执行重排
    reordered_df = reorder_hamming_matrix(csv_file, new_order)
    
    # 保存结果（可替换为其他文件名）
    reordered_df.to_csv(output_csv)
    print("矩阵重排完成，结果已保存为'reordered_mouse_hamming_matrix.csv'")
