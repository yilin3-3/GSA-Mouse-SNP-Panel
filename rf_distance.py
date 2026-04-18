from ete3 import Tree
import os
import numpy as np


def calculate_rf_distance_matrix(tree_files, output_file="rf_distance_matrix.txt"):
    """
    计算多个进化树之间的 Robinson-Foulds (RF) 距离矩阵。

    参数:
    tree_files: 进化树文件路径的列表，文件格式应为 Newick 格式。
    output_file: 输出距离矩阵的文件名，默认为 "rf_distance_matrix.txt"。

    返回:
    rf_matrix: 包含 RF 距离的二维列表。
    """

    # 检查输入文件是否存在
    for file in tree_files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"文件 {file} 不存在。")

    num_trees = len(tree_files)
    if num_trees < 2:
        raise ValueError("至少需要提供两个树文件来计算距离。")

    # 读取所有树并存储在列表中
    trees = []
    tree_names = []
    print("正在读取树文件...")
    for file in tree_files:
        try:
            # 使用 quoted_node_names=True 处理带特殊字符的节点名
            t = Tree(file, quoted_node_names=True)
            trees.append(t)
            # 使用文件名作为树的标识
            tree_names.append(os.path.basename(file))
            print(f"成功读取: {os.path.basename(file)}")
        except Exception as e:
            print(f"读取文件 {os.path.basename(file)} 失败: {e}")

    # 初始化距离矩阵
    rf_matrix = np.zeros((num_trees, num_trees))

    # 计算每对树之间的 RF 距离
    print("\n正在计算 RF 距离矩阵...")
    for i in range(num_trees):
        for j in range(i, num_trees):
            if i == j:
                rf_matrix[i, j] = 0
                rf_matrix[j, i] = 0
            else:
                # 计算两棵树之间的 RF 距离
                # 关键修改：修正括号匹配问题
                rf, max_rf, common_leaves, _, _, _, _ = trees[i].robinson_foulds(
                #rf, max_rf, common_leaves,parts_t1, parts_t2  = trees[i].robinson_foulds(
                trees[j], unrooted_trees=True
                )

                # 存储非标准化的 RF 距离
                rf_matrix[i, j] = rf
                rf_matrix[j, i] = rf  # 距离矩阵是对称的

                # 打印进度信息
                tree1_name = os.path.basename(tree_files[i])
                tree2_name = os.path.basename(tree_files[j])
                print(f"{tree1_name} vs {tree2_name}: RF = {rf}, 最大可能 RF = {max_rf}")

                # 将距离矩阵保存到文件
                with open(output_file, 'w') as f:
                # 写入表头（树文件名）
                    f.write("\t" + "\t".join(tree_names) + "\n")
                    for i in range(num_trees):
                # 写入每行的内容
                        row = [str(int(rf_matrix[i, j])) for j in range(num_trees)]
                        f.write(tree_names[i] + "\t" + "\t".join(row) + "\n")

                print(f"\nRF 距离矩阵已保存到: {output_file}")
    return rf_matrix


# -----------------------------------------------------------------------------
# 使用示例
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    # 1. 首先安装 ete3 库（如果还没安装）
    # !pip install ete3 numpy

    # 2. 准备你的进化树文件
    # 这里假设你有三个 Newick 格式的树文件
    # 请将下面的文件名替换成你自己的树文件路径

    tree_files = [
        #"mgp_REL2021_46_LDp4v7RA.nwk",
        "mgp_REL2021_LdW1p4v7R_L315669.nwk",
        "MGI_CS_phy_L512_11121025.nwk"  # 假设你也修复了这个文件
    ]

    # 3. 运行函数计算距离矩阵
    try:
        rf_matrix = calculate_rf_distance_matrix(tree_files)

        # 4. (可选) 直接打印矩阵查看结果
        print("\nRF 距离矩阵 (非标准化):")
        print(rf_matrix)

    except Exception as e:
        print(f"程序运行出错: {e}")
