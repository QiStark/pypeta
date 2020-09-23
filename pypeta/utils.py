import pandas as pd
import numpy as np
import os
import json


def is_float(string: str = '') -> bool:
    """Return True if the input string can be a float number

    Args:
        string (str, optional): input string. Defaults to ''.

    Returns:
        bool: return True if the input string can be converted to float 
    """
    try:
        float(string)
        return True
    except:
        return False


def restrict_series_value_counts_to_designated_records(ser: pd.Series,
                                                       limit: int = 20
                                                       ) -> pd.DataFrame:
    """Cut a long pandas Series to degisted records, the remain records will be aggregated in the 'Others' record.

    Args:
        ser (pd.Series): The long pandas Series to be processed
        limit (int, optional): Records number intended to be reamined. Defaults to 20.

    Returns:
        pd.DataFrame: processed records with the input Series`s index as the first column
    """
    length = len(ser)
    if length > limit:
        thres = limit - 1
        others = pd.Series(ser[thres:].sum(), index=['Others'])
        #ser = pd.concat([ser[:thres], others])
        ser = ser[:thres].append(others)

    df = pd.DataFrame(ser).reset_index()
    df.columns = pd.Index(['类别', '数量'])

    return df


def positive_rate(values: list, positive_tags: list) -> tuple:
    """Calculate positive tags marked values percentage in the total ones

    Args:
        values (list): the total values
        positive_tags (list): values that are regarded as positive values

    Returns:
        tuple: tuple for total values number, effective values number and percentage of positive values in the input values
    """
    values = list(values)

    total_value_num = len(values)
    missing_value_num = values.count(np.nan)
    effective_value_num = total_value_num - missing_value_num
    positvie_event_num = sum([values.count(tag) for tag in positive_tags])

    positive_rate = 0 if effective_value_num == 0 else positvie_event_num / effective_value_num

    return (total_value_num, effective_value_num, positive_rate)


def mut_freq_per_gene(maf_df: pd.DataFrame,
                      cnv_df: pd.DataFrame = pd.DataFrame([]),
                      sv_df: pd.DataFrame = pd.DataFrame([])) -> pd.Series:
    """Calculate variation frequency for each Gene in the input data set

    Args:
        maf_df (pd.DataFrame): maf records for SNV and InDel variations
        cnv_df (pd.DataFrame, optional): cnv records. Defaults to pd.DataFrame([]).
        sv_df (pd.DataFrame, optional): SV/Fusion/Gene rearrengement record. Defaults to pd.DataFrame([]).

    Raises:
        ValueError: input data set has no records

    Returns:
        pd.Series: variation freq for each Gene
    """
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


def filter_description(json_str: str) -> str:
    """Parse Peta restricts as literal description

    Args:
        json_str (str): string format Peta restricts

    Returns:
        str: literal description
    """
    filter_dict = json.loads(json_str)

    literal_description = ''
    literal_description += f'选取的研究数据集包括'
    literal_description += ','.join([*filter_dict['studyIds']])
    literal_description += '。'

    attributesRangeFilters = filter_dict['attributesRangeFilters']
    attributesEqualFilters = filter_dict['attributesEqualFilters']
    if attributesRangeFilters or attributesEqualFilters:
        literal_description += ','.join([*attributesRangeFilters,
                                        *attributesEqualFilters])

    return literal_description


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