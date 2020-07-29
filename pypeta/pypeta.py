import os
import sys
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import requests
import json
import time


class AccountError(RuntimeError):
    def __init__(self, account: str, message: str):
        self.account = account
        self.message = message


class NetworkError(RuntimeError):
    def __init__(self, message: str):
        self.message = message


class FetchError(RuntimeError):
    def __init__(self, restriction: str):
        self.restriction = restriction


class Peta(object):
    '''peta code interface class'''
    def __init__(self, username: str = '', password: str = '',
                 token: str = '', host: str='https://peta.bgi.com/api'):
        '''
        login
        '''
        self.host=host.strip('/')
        if token:
            jar = requests.cookies.RequestsCookieJar()
            jar.set('token', token)
            self.cookies = jar
        elif username and password:
            data = {'name': username, 'password': password}
            r = requests.post(f'{self.host}/peta/user/getticket',
                              data=data)
            if r.status_code != 200:
                raise NetworkError('login')
            elif r.text == '{}':
                raise AccountError(
                    username,
                    'Failed to access PETA. Error in username or password')
            else:
                self.cookies = r.cookies
        else:
            raise AccountError(
                ' Please specify token or username with password to access PETA'
            )

        self.headers = {
            "User-Agent":
            "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.132 Safari/537.36",
            'content-type': 'application/json'
        }

        self.data_restriction = {
            "studyIds": ["chol_nus_2012"],
            "attributesRangeFilters": [],
            "attributesEqualFilters": [],
            "mutationFilter": {
                "hugoGeneSymbols": [],
                "exacStart": 0,
                "exadEnd": 1,
                "vabundStart": 0,
                "vabundEnd": 1,
                "variantSource": [],
                "variantType": [],
                "variantClass": [],
                "sequencer": [],
                "sequencerSource": []
            },
            "cnvFilter": {},
            "svFilter": {},
            "pageIndex": 1,
            "pageSize": 999999,
        }

    def set_data_restriction_from_json_file(self, json_file: str):
        ''''''
        with open(json_file) as jsonin:
            data_restriction = json.load(jsonin)
            self.data_restriction = data_restriction

    def set_data_restriction_from_json_string(self, json_string: str):
        ''''''
        self.data_restriction = json.loads(json_string)

    # 4 feteh data interface, including clinical, mutation, cnv and sv. return dataframe
    def _fetch_data(self, url: str):
        '''fetch data from peta'''
        r = requests.post(url,
                          data=json.dumps(self.data_restriction),
                          cookies=self.cookies,
                          headers=self.headers)
        if r.status_code != 200:
            raise NetworkError('fetch data')
        elif re.findall(r'"responseCode":"-2"', r.text):
            raise FetchError(r.text)
        else:
            return r.text

    def fetch_mutation_data(self):
        url = f'{self.host}/peta/mutation/getMAFData'
        return pd.read_json(self._fetch_data(url))

    def fetch_clinical_data(self):
        url = f'{self.host}/peta/clinical/sampleClinicalData'

        #print("http start")
        #print(time.asctime())
        fetched_json = self._fetch_data(url)
        #print("http end\nread json start")
        #print(time.asctime())
        json_dict = json.loads(fetched_json)
        #print("read json end\ntable start")
        #print(time.asctime())

        samples_dict = json_dict['data']['samples']

        target_samples_list = []
        for sample_dict in samples_dict:

            target_sample_dict = {}
            for key, value in sample_dict.items():
                if key == 'clinicalData':
                    attr_dict_list = sample_dict['clinicalData']
                    for attr_dict in attr_dict_list:
                        for sub_key, sub_value in attr_dict.items():
                            target_sample_dict[
                                attr_dict['attrId']] = attr_dict['attrValue']
                else:
                    target_sample_dict[key] = value

            target_samples_list.append(target_sample_dict)

        df = pd.DataFrame(target_samples_list)

        #print("table end\n")
        #print(time.asctime())

        return df

    def fetch_cnv_data(self):
        return None
        url = f'{self.host}/peta/mutation/getCNVData'
        return self._fetch_data(url)

    def fetch_sv_data(self):
        return None
        url = f'{self.host}/peta/mutation/getSVData'
        return self._fetch_data(url)

    def fetch(self):
        '''fetch all info of the selected samples'''
        return [
            self.fetch_clinical_data(),
            self.fetch_mutation_data(), self.fetch_cnv_data, self.fetch_sv_data
        ]

    # list all the studys current user can see
    def list_visible_studys(self):
        r = requests.post(f'{self.host}/peta/home/getStudies',
                          data=json.dumps({
                              "name": "",
                              "parentType": [],
                              "groups": [],
                              "tags": []
                          }),
                          cookies=self.cookies,
                          headers=self.headers)
        if r.status_code != 200:
            raise NetworkError('list visible studys')
        elif re.findall(r'"responseCode":"-2"', r.text):  # 这里可能需要做单独的处理
            raise FetchError(self.data_restriction)
        else:
            study_list_json = json.loads(r.text)

            dfs = []
            for cancer_type in study_list_json['data']:

                for cancer_type_detail in study_list_json['data'][cancer_type][
                        'studies']:
                    df = pd.DataFrame(study_list_json['data'][cancer_type]
                                      ['studies'][cancer_type_detail]['data'])
                    df['cancerType'] = cancer_type
                    df['cancerTypeDetail'] = cancer_type_detail
                    dfs.append(df)

            return pd.concat(dfs, sort=True)

    # select studys
    def select_studys(self, study_ids: list = []):
        self.data_restriction['studyIds'] = study_ids

    # list all attributes of selected studys
    def list_sample_attributes(self):
        pass

    # designate attributes filters
    def designate_sample_filters(self, filters: defaultdict = {}):
        pass

    # list all variation atrributes of selected samples
    def list_variation_attributes(self):
        pass

    # set variation thresholds
    def set_variation_thresholds(self, thresholds: defaultdict = {}):
        pass

    # cnv和sv也增加

    # 样本 变异等多种查询
    def querys(self):
        pass

    # 多种分析 整合到一起

    # from json

    # from csv 增加一个文件接口

    def analysises(self):
        pass

    # 辅助函数
    def maf_to_yj(self):
        hzm = self.fetch_mutation_data()
        hzm = hzm[[
            'Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome',
            'Reference_Allele', 'Tumor_Seq_Allele2', 'Start_Position',
            'End_Position'
        ]]
        hzm = hzm.rename(
            columns={
                'Tumor_Sample_Barcode': 'SampleId',
                'Hugo_Symbol': 'Gene',
                'Chromosome': 'CHR',
                'Reference_Allele': 'REF',
                'Tumor_Seq_Allele2': 'ALT',
            })
        hzm['HGVSp'] = np.nan
        hzm['Case_Read1'] = 1
        hzm['Case_Read2'] = 1
        hzm['Case_Var_Freq'] = 1
        hzm = hzm.reindex(columns=[
            'SampleId', 'Gene', 'HGVSp', 'CHR', 'REF', 'ALT', 'Start_Position',
            'End_Position', 'Case_Read1', 'Case_Read2', 'Case_Var_Freq'
        ])

        return hzm
