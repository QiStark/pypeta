import os
import sys
import re
import numpy as np
import pandas as pd
from collections import defaultdict


class peta(object):
    '''Define peta class'''

    # only effective account can create the object
    def __init__(self,username:str,passcode:str):
        pass

    # login in
    def _login(self):
        pass
    # return token?

    # get full data for selected samples, 4 dataframes including clinical, mutation, cnv and sv 
    def fetch_data(self):
        pass

    # list all the studys current user can see
    def list_visible_studys(self):
        pass

    # select studys 
    def select_studys(self,study_names:list=[]):
        pass

    # list all attributes of selected studys
    def list_sample_attributes(self):
        pass

    # designate attributes filters
    def designate_sample_filters(self,filters:defaultdict={}):
        pass

    # list all variation atrributes of selected samples
    def list_variation_attributes(self):
        pass

    # set variation thresholds
    def set_variation_thresholds(self, thresholds:defaultdict={}):
        pass
    
    # 样本 变异等多种查询
    def querys(self):
        pass


    # 多种分析 整合到一起
    def analysises(self):
        pass