import os
import sys
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import requests
import json


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
    def __init__(self, username: str, password: str):
        '''
        login
        '''
        data = {'name': username, 'password': password}
        r = requests.post('https://peta.bgi.com/api/peta/user/getticket',
                          data=data)
        if r.status_code != 200:
            raise NetworkError('login')
        elif r.text == '{}':
            raise AccountError(username, 'Error in username or password')
        else:
            self.cookies = r.cookies

        self.headers = {
            "User-Agent":
            "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.132 Safari/537.36",
            'content-type': 'application/json'
        }

        self.data_restriction = {
            "studyIds": ["chol_nus_2012"],
            "attributesRangeFilters": [],
            "attributesEqualFilters": [{
                "attributeId": "OS_STATUS",
                "attributeType": "PATIENT",
                "values": ["ALIVE"]
            }],
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
            "svFilter": {}
        }

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
            raise FetchError
        else:
            return r.text

    def fetch_mutation_data(self):
        url = 'https://peta.bgi.com/api/peta/mutation/getMAFData'
        return pd.read_json(self._fetch_data(url))

    def fetch_clinical_data(self):
        url = 'https://peta.bgi.com/api/peta/clinical/sampleClinicalData'

        samples_data = pd.read_json(self._fetch_data(url)).iloc[1, 0]

        dfs = []
        for sample_record in samples_data:
            sampleId = sample_record['sampleId']
            patientId = sample_record['patientId']

            df = pd.DataFrame(
                sample_record['clinicalData']).set_index('attrId').T.iloc[0:1]
            df['patientId'] = patientId
            df['sampleId'] = sampleId

            dfs.append(df)

        return pd.concat(dfs)

    def fetch_cnv_data(self):
        url = 'https://peta.bgi.com/api/peta/mutation/getCNVData'
        return self._fetch_data(url)

    def fetch_sv_data(self):
        url = 'https://peta.bgi.com/api/peta/mutation/getSVData'
        return self._fetch_data(url)

    # list all the studys current user can see
    def list_visible_studys(self):
        r = requests.post('https://peta.bgi.com/api/peta/home/getStudies',
                          data=json.dumps({"name":"","parentType":[],"groups":[],"tags":[]}),
                          cookies=self.cookies,
                          headers=self.headers)
        if r.status_code != 200:
            raise NetworkError('list visible studys')
        elif re.findall(r'"responseCode":"-2"', r.text):  # 这里可能需要做单独的处理
            raise FetchError(self.data_restriction)
        else:
            study_list_json=json.loads(r.text)

            dfs=[]
            for cancer_type in study_list_json['data']:

                for cancer_type_detail in study_list_json['data'][cancer_type]['studies']:
                    df=pd.DataFrame(study_list_json['data'][cancer_type]['studies'][cancer_type_detail]['data'])
                    df['cancerType']=cancer_type
                    df['cancerTypeDetail']=cancer_type_detail
                    dfs.append(df)

            return pd.concat(dfs,sort=True)


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
