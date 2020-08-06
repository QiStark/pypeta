import pandas as pd
import numpy as np


def is_float(string: str = ''):
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
        ser = pd.concat([ser[:thres], others])

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
