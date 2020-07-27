# Pypeta
pypeta is the python code interface for BGI-PETA database. It provides APIs for login, querying and fetching data.

## Installation
From source code.
```
$ git clone https://github.com/JaylanLiu/pypeta
$ cd pypeta
$ python setup.py install 
```

## Example usage
Import libraies:
```
import pypeta
```

Login to BGI-PETA database via PETA account:
```
peta=pypeta.Peta(username='bgi-peta@genomics.cn',password='******')
```
Or login via token:
```
peta=pypeta.Peta(token = '******')
```

Get all visable studies info:
```
peta.list_visible_studys()
```

Set target study id and other filters:
```
peta.select_studys(["chol_nus_2012"])
peta.data_restriction["attributesEqualFilters"]=[]
```

Also studies and filters can be load from a json file, which can be get from the PETA website:
```
peta.set_data_restriction_from_json_file('json.txt')
```

Fetch SNV and InDel data as MAF format:
```
peta.fetch_mutation_data()
```

Fetch clinical info as tsv format:
```
peta.fetch_clinical_data()
```