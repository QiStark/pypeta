import pandas as pd
import numpy as np
import os
import json


def is_float(string: str = '') -> bool:
    '''return True if a string can be converted to float'''
    try:
        float(string)
        return True
    except:
        return False


def restrict_series_value_counts_to_designated_records(ser: pd.Series,
                                                       limit: int = 20):
    '''
    组合的index用什么和列标签的命名需要考虑中英文
    '''
    length = len(ser)
    if length > limit:
        thres = limit - 1
        others = pd.Series(ser[thres:].sum(), index=['Others'])
        #ser = pd.concat([ser[:thres], others])
        ser = ser[:thres].append(others)

    df = pd.DataFrame(ser).reset_index()
    df.columns = pd.Index(['类别', '数量'])

    return df


def positive_rate(values: list, positive_tags: list):
    values = list(values)

    total_value_num = len(values)
    missing_value_num = values.count(np.nan)
    effective_value_num = total_value_num - missing_value_num
    positvie_event_num = sum([values.count(tag) for tag in positive_tags])

    positive_rate = 0 if effective_value_num == 0 else positvie_event_num / effective_value_num

    return (total_value_num, effective_value_num, positive_rate)


def mut_freq_per_gene(maf_df: pd.DataFrame,
                      cnv_df: pd.DataFrame = pd.DataFrame([]),
                      sv_df: pd.DataFrame = pd.DataFrame([])):
    '''使用计算每个基因在群体中的突变频率'''
    mut_df = maf_df[['Tumor_Sample_Barcode', 'Hugo_Symbol']]
    mut_df.columns = pd.Index(['Sample_ID', 'Hugo_Symbol'])

    if len(cnv_df) == 0:
        pass

    if len(sv_df) == 0:
        pass

    samples_num = len(mut_df.Sample_ID.drop_duplicates())
    if samples_num == 0:
        raise ValueError

    return mut_df.dropna().drop_duplicates().Hugo_Symbol.value_counts(
    ) / samples_num


def filter_description(json_str: str):
    """Parse Peta restricts as literal description

    Args:
        json_str (str): string format Peta restricts
    """
    filter_dict = json.loads(json_str)
    print(f'选取的研究数据集包括', end='')
    print(*filter_dict['studyIds'], sep=',', end='')
    print('。')

    attributesRangeFilters = filter_dict['attributesRangeFilters']
    attributesEqualFilters = filter_dict['attributesEqualFilters']
    if attributesRangeFilters or attributesEqualFilters:
        print('样本过滤条件为',
              *attributesRangeFilters,
              *attributesEqualFilters,
              sep=',')


def get_certain_file_type_from_certain_depth_folders(root_dir: str,
                                                     suffix: list,
                                                     depth: int = 1) -> list:
    """Return list of file path with specified suffix in the specified depth of input directory

    Args:
        root_dir (str): input derectory
        suffix (list): target file suffix, like xlsx
        depth (int, optional): sub-folder depth to search. Defaults to 1.

    Returns:
        list: target file paths list
    """

    branch_paths = [root_dir]
    for d in range(depth - 1):
        leaf_paths = branch_paths
        branch_paths = []
        for leaf_path in leaf_paths:
            branch_paths.extend([
                os.path.join(leaf_path, sub_f)
                for sub_f in os.listdir(leaf_path)
                if os.path.isdir(os.path.join(leaf_path, sub_f))
            ])
    else:
        leaf_paths = branch_paths

    target_files = []
    for leaf_path in leaf_paths:
        target_files.extend([
            os.path.join(leaf_path, sub_f) for sub_f in os.listdir(leaf_path)
            if sub_f.split('.')[-1] in suffix
        ])

    return target_files