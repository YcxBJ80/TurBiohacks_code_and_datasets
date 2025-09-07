#!/usr/bin/env python3
"""
数据集合并脚本
合并 merged_features_labels.csv 和 GSE25097.csv
只保留两个数据集的交集特征（列），并调整标签编码
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from sklearn.model_selection import train_test_split

def merge_datasets(file1_path, file2_path, output_path):
    """
    合并两个CSV数据集，只保留交集特征
    
    Args:
        file1_path: merged_features_labels.csv 路径
        file2_path: GSE25097.csv 路径  
        output_path: 合并后的输出文件路径
    """
    
    print("🔄 开始合并数据集...")
    
    # 读取第一个数据集 (merged_features_labels.csv)
    print(f"📊 读取第一个数据集: {file1_path}")
    df1 = pd.read_csv(file1_path, index_col=0)
    print(f"   - 样本数: {len(df1)}")
    print(f"   - 特征数: {len(df1.columns) - 1} (不含Label)")
    
    # 读取第二个数据集 (GSE25097.csv)
    print(f"📊 读取第二个数据集: {file2_path}")
    df2 = pd.read_csv(file2_path, index_col=0)
    print(f"   - 样本数: {len(df2)}")
    print(f"   - 特征数: {len(df2.columns) - 1} (不含Label)")
    
    # 检查标签分布
    print("\n🏷️ 原始标签分布:")
    print("数据集1 (merged_features_labels.csv):")
    if 'Label' in df1.columns:
        label_col1 = 'Label'
    else:
        # 假设最后一列是标签列
        label_col1 = df1.columns[-1]
    print(f"   {df1[label_col1].value_counts().to_dict()}")
    
    print("数据集2 (GSE25097.csv):")
    if 'Label' in df2.columns:
        label_col2 = 'Label'
    else:
        # 假设最后一列是标签列
        label_col2 = df2.columns[-1]
    print(f"   {df2[label_col2].value_counts().to_dict()}")
    
    # 分离特征和标签
    features1 = df1.drop(columns=[label_col1])
    labels1 = df1[label_col1]
    
    features2 = df2.drop(columns=[label_col2])  
    labels2 = df2[label_col2]
    
    # 找到交集特征
    common_features = list(set(features1.columns) & set(features2.columns))
    print(f"\n🔗 找到交集特征: {len(common_features)} 个")
    print(f"   数据集1独有特征: {len(features1.columns) - len(common_features)} 个")
    print(f"   数据集2独有特征: {len(features2.columns) - len(common_features)} 个")
    
    if len(common_features) == 0:
        raise ValueError("❌ 错误：两个数据集没有共同的特征列！")
    
    # 只保留交集特征
    features1_common = features1[common_features]
    features2_common = features2[common_features]
    
    # 处理标签编码
    print("\n🔄 处理标签编码...")
    
    # 数据集1: 将字符串标签转换为数值
    # "Normal" -> 0, "HBV-HCC" -> 1
    labels1_numeric = labels1.map({
        "Normal": 0,
        "HBV-HCC": 1
    })
    
    # 数据集2: 将标签1转换为2，保持0不变
    # 0 -> 0, 1 -> 2  
    labels2_numeric = labels2.map({0: 0, 1: 2})
    
    print("处理后标签分布:")
    print(f"   数据集1: {dict(labels1_numeric.value_counts().sort_index())}")
    print(f"   数据集2: {dict(labels2_numeric.value_counts().sort_index())}")
    
    # 重新添加标签列
    features1_common['Label'] = labels1_numeric
    features2_common['Label'] = labels2_numeric
    
    # 合并数据集
    print(f"\n🔀 合并数据集...")
    merged_df = pd.concat([features1_common, features2_common], axis=0, ignore_index=False)
    
    print(f"✅ 合并完成:")
    print(f"   - 总样本数: {len(merged_df)}")
    print(f"   - 特征数: {len(merged_df.columns) - 1} (不含Label)")
    print(f"   - 最终标签分布: {dict(merged_df['Label'].value_counts().sort_index())}")
    
    # 保存合并后的数据集
    print(f"\n💾 保存完整数据集到: {output_path}")
    merged_df.to_csv(output_path, index=True)
    
    # 分割训练集和验证集 (8:2)
    print("\n🔄 按8:2比例分割训练集和验证集...")
    
    # 为了保持类别平衡，使用分层抽样
    X = merged_df.drop('Label', axis=1)
    y = merged_df['Label']
    
    # 分层抽样确保训练集和验证集中各类别比例相同
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    # 重新组合特征和标签
    train_df = X_train.copy()
    train_df['Label'] = y_train
    
    val_df = X_val.copy()
    val_df['Label'] = y_val
    
    # 保存训练集和验证集
    train_path = str(output_path).replace('.csv', '_train.csv')
    val_path = str(output_path).replace('.csv', '_valid.csv')
    
    print(f"💾 保存训练集到: {train_path}")
    train_df.to_csv(train_path, index=True)
    
    print(f"💾 保存验证集到: {val_path}")
    val_df.to_csv(val_path, index=True)
    
    # 生成统计报告
    report = f"""
=== 数据集合并和分割报告 ===
输入文件:
  - 数据集1: {file1_path}
  - 数据集2: {file2_path}

合并结果:
  - 完整数据集: {output_path} ({len(merged_df)} 样本)
  - 训练集: {train_path} ({len(train_df)} 样本, {len(train_df)/len(merged_df)*100:.1f}%)
  - 验证集: {val_path} ({len(val_df)} 样本, {len(val_df)/len(merged_df)*100:.1f}%)
  - 共同特征数: {len(common_features)}
  
标签编码:
  - 0: 正常/非肿瘤样本
  - 1: HBV-HCC样本  
  - 2: 其他肿瘤样本

标签分布:
  完整数据集: {dict(merged_df['Label'].value_counts().sort_index())}
  训练集: {dict(train_df['Label'].value_counts().sort_index())}
  验证集: {dict(val_df['Label'].value_counts().sort_index())}

数据来源:
  - Label 0: 来自两个数据集的正常/非肿瘤样本
  - Label 1: 来自 merged_features_labels.csv 的 HBV-HCC 样本
  - Label 2: 来自 GSE25097.csv 的肿瘤样本
"""
    
    print(report)
    
    # 保存报告
    report_path = str(output_path).replace('.csv', '_report.txt')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    
    return merged_df, train_df, val_df

def main():
    """主函数"""
    # 设置文件路径
    current_dir = Path(__file__).parent
    file1_path = current_dir / "merged_features_labels.csv"
    file2_path = current_dir / "GSE25097.csv"
    output_path = current_dir / "combined_dataset.csv"
    
    # 检查输入文件是否存在
    if not file1_path.exists():
        raise FileNotFoundError(f"❌ 文件不存在: {file1_path}")
    if not file2_path.exists():
        raise FileNotFoundError(f"❌ 文件不存在: {file2_path}")
    
    try:
        # 执行合并和分割
        merged_df, train_df, val_df = merge_datasets(file1_path, file2_path, output_path)
        print("\n🎉 数据集合并和分割成功完成！")
        
    except Exception as e:
        print(f"❌ 合并过程中出错: {str(e)}")
        raise

if __name__ == "__main__":
    main()
