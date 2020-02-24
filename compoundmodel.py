# import numpy as np
import base64
import pickle

import pandas as pd
from time import time
# import matplotlib.pyplot as plt
import os
# import shutil
import ast
import math

import logging
import requests
from werkzeug.utils import secure_filename

logging.getLogger("requests").setLevel(logging.ERROR)
import xml.etree.ElementTree as ET
from tqdm import tqdm
from flask import Flask
import plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
import json
import glob
# import seaborn as sns
# from pandas.plotting import scatter_matrix
from scipy import stats
from plotly.subplots import make_subplots

#from IPython.display import Image
# from sklearn.metrics import f1_score
# from sklearn.metrics import confusion_matrix
# import seaborn as sns
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import MolFromSmiles, MolToSmiles
from padelpy import padeldescriptor
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
# from functools import reduce
# from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
# from sklearn.metrics import matthews_corrcoef, precision_recall_fscore_support, roc_curve, auc, confusion_matrix
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn import model_selection, metrics
# from sklearn import model_selection, tree, preprocessing, metrics, linear_model
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.naive_bayes import GaussianNB
# from sklearn.tree import DecisionTreeClassifier
# from sklearn.ensemble import GradientBoostingClassifier
from sklearn.feature_selection import RFECV
# from sklearn.feature_selection import RFE, RFECV
from sklearn.linear_model import LogisticRegression
import time
# from hpsklearn import HyperoptEstimator, one_vs_rest, random_forest
# import hpsklearn.demo_support
# from hyperopt import tpe
import numpy as np
# from sklearn.metrics import roc_auc_score
# import matplotlib.pyplot as plt
# from sklearn.datasets import make_classification
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedShuffleSplit
from imblearn.combine import SMOTEENN
import shutil
import re
from flask import request
from biothings_client import get_client

# Application
dir = os.path.abspath(os.getcwd())
dir = (dir+'\\file\\')
dir = dir.replace('\\', '/')
server = Flask(__name__)
server.config['UPLOAD_FOLDER'] = dir


def gettargetdetail(uniID):
    detail = {'Function': '-', 'Subunit': '-', 'Similarity': '-', 'URL': '-'}
    urluniprot = '{http://uniprot.org/uniprot}'
    uniprotxml = requests.get("https://www.uniprot.org/uniprot/"+uniID+".xml")
    if uniprotxml.text != '':
        root = ET.fromstring(uniprotxml.text)
        rv = root[0].attrib['dataset']
        sim = root.find(urluniprot + 'entry')
        if sim is not None:
            for child in sim:
                if child.tag == urluniprot + 'comment':
                    if child.attrib['type'] == 'function':
                        detail['Function'] = child.find(urluniprot + 'text').text
                        print(child.find(urluniprot + 'text').text)
                    if child.attrib['type'] == 'subunit':
                        detail['Subunit'] = child.find(urluniprot + 'text').text
                        print(child.find(urluniprot + 'text').text)
                    if child.attrib['type'] == 'similarity':
                        detail['Similarity'] = child.find(urluniprot + 'text').text
                        print(child.find(urluniprot + 'text').text)
    detail['URL'] = 'http://uniprot.org/uniprot/'+uniID
    return detail


# Convert IC50 to pIC50
from math import log10
def pIC50(df0):
    pIC50 = []
    for i in df0.value:
        # Converts nM to M
        molar = i * (10 ** -9)
        pIC50.append(-log10(molar))
    df0['pIC50'] = pIC50
    df0 = df0.drop('value', 1)
    return df0


# Define active, inactive, intermediate compounds
def STATUS(df):
    STATUS = []
    count = 0
    avg = 0
    for i in df.pIC50:
        avg = avg + i
        count = count + 1
    avg = avg / count
    print("Cut-off = " + str(avg))
    for i in df.pIC50:
        if i > avg:
            STATUS.append("active")  # active
        # elif 1000 < i < 10000:
        # STATUS.append("intermediate") #intermediate
        elif i <= avg:
            STATUS.append("inactive")  # inactive
    df['STATUS'] = STATUS
    return df


def select_status(df):
    from chembl_webresource_client.new_client import new_client
    #     print("##########################################################")
    #     print('For the ' + str(df.name) + ' dataset...')
    #     print('')
    active1 = df[(df['STATUS'] == 'active') &
                 (df['relation'] == '=')]
    active2 = df.loc[df.index.isin(
        df[(df['relation'] == '<') &
           (df['STATUS'] == 'active')].index)]
    Active = active1.append(active2)
    Active_dup = Active[Active['molecule_chembl_id'].duplicated()]
    #     print ("For total of : "+ str(len(df)) + ' compounds')
    #     print('')
    #     print("Active dataset: "+ str(len(Active)) + ' with ' + str(len(Active_dup)) + ' duplicate compounds')
    inactive1 = df[(df['STATUS'] == 'inactive') &
                   (df['relation'] == '=')]
    inactive2 = df.loc[df.index.isin(
        df[(df['relation'] == '>') &
           (df['STATUS'] == 'inactive')].index)]
    Inactive = inactive1.append(inactive2)
    Inactive_dup = Inactive[Inactive['molecule_chembl_id'].duplicated()]
    #     print("Inactive dataset: "+ str(len(Inactive))+ ' with ' + str(len(Inactive_dup)) + ' duplicate compounds')
    Active_can_use = Active.loc[~Active['molecule_chembl_id'].isin(Inactive['molecule_chembl_id'])]
    Inactive_can_use = Inactive.loc[~Inactive['molecule_chembl_id'].isin(Active['molecule_chembl_id'])]
    #     print("There are: "+ str(len(Active_can_use)) + ' compounds can be used in ' + str(len(Active))+ ' Active dataset')
    #     print("There are: "+ str(len(Inactive_can_use)) + ' compounds can be used in ' + str(len(Inactive))+ ' Inactive dataset')
    #     print(str(len(Active_can_use)+len(Inactive_can_use)))
    df2 = Active_can_use.append(Inactive_can_use)
    #     print(str(len(df2)) + ' compounds can be uesd' + ' from ' + str(len(df)))
    return df2


def Chembl_Activity_REST(chemblID):
    actdf = pd.DataFrame()
    resp = requests.get('https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id__exact=' + chemblID + '&limit=1000')
    if resp.status_code != 200:
        # This means something went wrong.
        print("Error: " + str((resp.status_code)))
    res = resp.json()
    actdf = pd.DataFrame(res.get('activities'))
    meta = res.get('page_meta')
    if meta['next']:
        url = ("https://www.ebi.ac.uk" + str(meta['next']))
    while(meta['next']):
        resp = requests.get(url)
        if resp.status_code != 200:
            # This means something went wrong.
            print("Error: " + str((resp.status_code)))
        res = resp.json()
        temp = pd.DataFrame(res.get('activities'))
        actdf = actdf.append(temp)
        meta = res.get('page_meta')
        url = ("https://www.ebi.ac.uk" + str(meta['next']))
    return actdf


def Chembl_Assay_REST(chemblID):
    actdf = pd.DataFrame()
    resp = requests.get('https://www.ebi.ac.uk/chembl/api/data/assay.json?target_chembl_id__exact=' + chemblID + '&limit=1000')
    if resp.status_code != 200:
        # This means something went wrong.
        print("Error: " + str((resp.status_code)))
    res = resp.json()
    actdf = pd.DataFrame(res.get('assays'))
    meta = res.get('page_meta')
    if meta['next']:
        url = ("https://www.ebi.ac.uk" + str(meta['next']))
    while(meta['next']):
        resp = requests.get(url)
        if resp.status_code != 200:
            # This means something went wrong.
            print("Error: " + str((resp.status_code)))
        res = resp.json()
        temp = pd.DataFrame(res.get('activities'))
        actdf = actdf.append(temp)
        meta = res.get('page_meta')
        url = ("https://www.ebi.ac.uk" + str(meta['next']))
    return actdf

def Load_compounds(df,IC50_ready):
    
    from chembl_webresource_client.new_client import new_client
    
    compound = pd.DataFrame()
    count = 0
    molecule = new_client.molecule
    keys = []
    bad_key = []
    for i in tqdm(df['molecule_chembl_id']):
        keys.append(i)
        count = count + 1
        if (count % 250) == 0 or count == len(df['molecule_chembl_id']):
            records = molecule.get(keys)
            keys = []
            temp = pd.DataFrame(records)
            compound = compound.append(temp)
    # for i in df['molecule_chembl_id']:
    #     try:
    #         c = new_client.molecule.get(i)
    #          print(type(c))
    #          c = pd.DataFrame(c)
    #         compound = compound.append(c,ignore_index=True)
    #     except:
    #                             bad_key.append(i)
    # print("Bad IDs: " + str(bad_key))
    # try:
    #     compound = new_client.molecule.get(list(df['molecule_chembl_id'].head(5)))
    # except:
    #                                        print(i['molecule_chembl_id'])
    # compound = new_client.molecule.get(list(df['molecule_chembl_id'].head(5)))
    # compound = new_client.molecule.get(['CHEMBL6498', 'CHEMBL6499', 'CHEMBL6505'])

    print("##########################################################")
    print('For the ' + str(df.name) + ' dataset...')
    print('')
    print('There are ' + str(len(df)) + ' bioactivity...')

    df2 = pd.DataFrame(compound)
    df3 = pd.merge(df2, IC50_ready, how='inner', on='molecule_chembl_id', suffixes=('', '_y'))
    df4 = df3.drop_duplicates(subset='molecule_chembl_id') # drop duplicate in chemblId
    #df5 = df3.drop_duplicates(subset='smiles')   # drop duplicate SMILES
    df6 = df4[pd.notnull(df4['canonical_smiles'])]         # drop Nan in SMILES
    
    print('get '+ str(len(df2)) + ' compounds')
    print('After remove duplicate ChEMBL ID there are  '+ str(len(df6)) + ' compounds')
    #print ('After remove duplicate and absent SMILES there are  '+ str(len(df6)) + ' compounds')
    print('')
    
    return df6


# clean smiles
def clean_smiles(df):
    remover = SaltRemover()
    SMILES_desalt = []
    #     print('Processing ' + str(df.name) + ' dataset...')
    for i in df.canonical_smiles:
        mol = MolFromSmiles(i)
        mol_desalt = remover.StripMol(mol)
        mol_SMILES = MolToSmiles(mol_desalt)
        SMILES_desalt.append(mol_SMILES)
    df['SMILES_desalt'] = SMILES_desalt
    df['mol'] = mol
    df2 = df
    return df2

def drop_duplicate (df):
    df2 = df.drop_duplicates(subset='SMILES_desalt')
    # print ("RAW data of " + str(len(df)) + \
    #       " SMILES has been reduced to "+ str(len(df2)) + " SMILES.")
    return df2


import matplotlib.pyplot as plt
import numpy


def plot_status(IC50, INH, path):
    # Calculate number of compounds
    n_compoundsIC50 = len(IC50)
    n_compoundsINH = len(INH)

    # Calculate number of features
    # n_features = All_data_ER_alpha.shape[1] -1

    # Calculate active compounds
    n_ActIC50 = len(IC50[IC50.STATUS == 'active'])
    n_ActINH = len(INH[INH.STATUS == 'active'])

    # Calculate inactive compounds
    n_InactIC50 = len(IC50[IC50.STATUS == 'inactive'])
    n_InactINH = len(INH[INH.STATUS == 'inactive'])

    # Calculate active compounds
    n_interIC50 = len(IC50[IC50.STATUS == 'intermediate'])
    n_interINH = len(INH[INH.STATUS == 'intermediate'])

    # Calculate graduation rate
    Active_rateIC50 = (n_ActIC50 / float(n_compoundsIC50)) * 100.
    Active_rateINH = (n_ActINH / float(n_compoundsINH)) * 100.
    # Print the results
    #     print ("Total number of compounds: {}".format(n_compounds))
    #     #print "Number of features: {}".format(n_features))
    #     print ("Number of Active compounds: {}".format(n_Act))
    #     print ("Number of Inactive compounds: {}".format(n_Inact))
    #     print ("Number of intermediate compounds: {}".format(n_inter))
    #     print ("Active compounds rate of the data: {:.2f}%".format(Active_rate))

    #     print('*******************************************')
    # ax = sns.countplot(x='STATUS', data=df)
    # fig = ax.get_figure()
    # fig.tight_layout()
    # fig.savefig(os.path.join(path, str(df.name) + '_Bar.pdf'), dpi=300)

    axisx = ['Active', 'Inactive']
    fig = go.Figure(data=[
        go.Bar(name='IC<sub>50</sub>', x=axisx, y=[n_ActIC50, n_InactIC50], marker_color='#ef6c00'),
        go.Bar(name='INH', x=axisx, y=[n_ActINH, n_InactINH], marker_color='#4db6ac')
    ])
    # Change the bar mode
    fig.update_layout(
        autosize=True,
        barmode='group',
        height=300,
        margin=dict(l=10, r=20, t=30, b=10)
    )
    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON


def getBioActs(filename, uniID):

    response = ""
    allgraph = []

    # if os.path.exists(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)):
    #     response = 'ok'
    #     return response

    pd.options.mode.chained_assignment = None

    path = os.path.join(server.config['UPLOAD_FOLDER'], filename+"_"+uniID)

    print("Working on: " + uniID)
    target = new_client.target.filter(target_components__accession=uniID)

    # Load targets in data frame
    targetsDF = pd.DataFrame(target)
    # Target doesn't exist in Chembl
    if targetsDF.empty:
        print("Bad Target: " + uniID)
        response = "Target doesn't not exist in Chembl."
        print(response)
        return response
    targetID = targetsDF.at[0, 'target_chembl_id']
    # print(targetID)

    # Get Bioactivities with targetID
    # bioactsDF = new_client.activity.filter(target_chembl_id=targetID)
    bioactsDF = Chembl_Activity_REST(targetID)
    # Get Assay with targetID, for confidence score
    # assayDF = new_client.assay.filter(target_chembl_id=targetID)
    assayDF = Chembl_Assay_REST(targetID)
    # Convert to DF
    bioactsDF = pd.DataFrame(bioactsDF)
    assayDF = pd.DataFrame(assayDF)
    # Bad Targets
    if len(bioactsDF) == 0 or len(assayDF) == 0:
        print("Bad Target: " + uniID)
        response = "Target has no associated bioactivity or binding type assay."
        print(response)
        return response

    # Merge both DFs
    bioactsDF = pd.merge(bioactsDF, assayDF[['assay_chembl_id', 'confidence_score']], how='inner', on='assay_chembl_id')

    # Drop values contain Empty values
    bioactsDF['value'].replace('', np.nan, inplace=True)
    bioactsDF.dropna(subset=['value'], inplace=True)
    # Drop operator contain Empty values
    bioactsDF['relation'].replace('', np.nan, inplace=True)
    bioactsDF.dropna(subset=['relation'], inplace=True)


    bioact = []

    labels = list(bioactsDF['confidence_score'].value_counts().index)
    values = list(bioactsDF['confidence_score'].value_counts().values)
    pulling = [0] * len(values)
    if 9 in labels:
        pulling[labels.index(9)] = 0.2
    # pull is given as a fraction of the pie radius
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, pull=pulling)])
    fig.update_layout(height=300, margin=dict(l=10, r=20, t=30, b=10), title_text='Confidence Score')

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    bioact.append(graphJSON)

    labels = list(bioactsDF['assay_type'].value_counts().index)
    values = list(bioactsDF['assay_type'].value_counts().values)
    pulling = [0] * len(values)
    if 'B' in labels:
        pulling[labels.index('B')] = 0.2
    # pull is given as a fraction of the pie radius
    fig2 = go.Figure(data=[go.Pie(labels=labels, values=values, pull=pulling)])
    fig2.update_layout(height=300, margin=dict(l=10, r=20, t=30, b=10), title_text='Assay Type')

    graphJSON = json.dumps(fig2, cls=plotly.utils.PlotlyJSONEncoder)
    bioact.append(graphJSON)

    labels = list(bioactsDF['type'].value_counts().index)
    values = list(bioactsDF['type'].value_counts().values)
    pulling = [0] * len(values)
    if 'IC50' in labels:
        pulling[labels.index('IC50')] = 0.2
    if 'INH' in labels:
        pulling[labels.index('INH')] = 0.2
    # pull is given as a fraction of the pie radius
    fig3 = go.Figure(data=[go.Pie(labels=labels, values=values, pull=pulling)])
    fig3.update_layout(height=300, margin=dict(l=10, r=20, t=30, b=10), title_text='Bioactivity Type')

    graphJSON = json.dumps(fig3, cls=plotly.utils.PlotlyJSONEncoder)
    bioact.append(graphJSON)

    # Making Directory
    if not os.path.exists(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)):
        os.makedirs(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID))

    with open(os.path.join(path, "_bioact.txt"), "w") as f:
        for s in bioact:
            f.write(str(s) + "\n")

    bioactsDF = bioactsDF[(bioactsDF['confidence_score'] == 9)]  # Direct single protein target assigned
    bioactsDF = bioactsDF[(bioactsDF['assay_type'] == 'B')]  # only binding data
    bioactsDF[['value']] = bioactsDF[['value']].apply(pd.to_numeric)

    # Prepare Dataset
    IC50 = bioactsDF[(bioactsDF['type'] == 'IC50')]  # 50% inhibition concentration
    Inhibition = bioactsDF[(bioactsDF['type'] == 'INH')]

    # Bad Targets
    if len(IC50) == 0 or len(Inhibition) == 0:
        print("Bad Target: " + uniID)
        response = "Target has no bioactivity type IC50 or INH."
        print(response)
        return response

    # Convert to pIC50
    IC50 = pIC50(IC50)
    IC50.name = 'IC50'
    Inhibition.name = 'INH'

    binmax = math.ceil(max(list(IC50['pIC50'].value_counts(sort=False).index)))
    label = []
    for k in range(0, binmax-1):
        label.append('('+str(k)+','+str(k+1)+']')
    IC50['bin'] = pd.cut(IC50['pIC50'], list(range(0, binmax)), labels=label)


    # Show Bar Chart before preprocessing
    # fig, ax = plt.subplots()
    # IC50['bin'].value_counts(sort=False).plot(ax=ax, kind='bar')
    # plt.xlabel('pIC50', size=12, fontweight='bold');
    # plt.ylabel('Frequency', size=12, fontweight='bold');
    # plt.tight_layout()
    # plt.savefig(os.path.join(path,'IC50_pIC50_Bar_Before.pdf'), dpi=300)

    indexx = list(IC50['bin'].value_counts(sort=False).index)
    valuesy = list(IC50['bin'].value_counts(sort=False).values)
    fig = go.Figure([
        go.Bar(
            x=indexx,
            y=valuesy,
            marker_color='#ef6c00'
        )
    ])
    fig.update_layout(
        autosize=True,
        height=300,
        margin=dict(l=10, r=20, t=30, b=10)
    )
    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    allgraph.append(graphJSON)

    IC50 = STATUS(IC50)
    IC50_ready = select_status(IC50)

    # Show Bar Chart after preprocessing
    # fig, ax = plt.subplots()
    # IC50_ready['bin'].value_counts(sort=False).plot(ax=ax, kind='bar')
    # plt.xlabel('pIC50', size=12, fontweight='bold');
    # plt.ylabel('Frequency', size=12, fontweight='bold');
    # plt.tight_layout()
    # plt.savefig(os.path.join(path, 'dataset pIC50 final.pdf'), dpi=300)

    indexx = list(IC50_ready['bin'].value_counts(sort=False).index)
    valuesy = list(IC50_ready['bin'].value_counts(sort=False).values)
    fig = go.Figure([
        go.Bar(
            x=indexx,
            y=valuesy,
            marker_color='#ef6c00'
        )
    ])
    fig.update_layout(
        autosize=True,
        height=300,
        margin=dict(l=10, r=20, t=30, b=10)
    )
    graphJSON2 = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    allgraph.append(graphJSON2)

    IC50 = STATUS(IC50)
    IC50_ready = select_status(IC50)

    # Preprecessing
    IC50_ready.name = 'IC50'
    IC50_raw = Load_compounds(IC50_ready, IC50_ready)
    Inhibition_raw = Load_compounds(Inhibition, IC50_ready)
    # Bad Targets
    if len(IC50_raw) == 0 or len(Inhibition_raw) == 0:
        # removing directory
        shutil.rmtree(path)
        print("Bad Target: " + uniID)
        response = "Target has no compound associated to activities."
        print(response)
        return response

    IC50_raw.name = 'I50_Raw'
    Inhibition_raw.name = 'INH_Raw'
    IC50_clean = clean_smiles(IC50_raw)
    Inhibition_clean = clean_smiles(Inhibition_raw)
    IC50_clean = drop_duplicate(IC50_clean)
    Inhibition_clean = drop_duplicate(Inhibition_clean)
    IC50_clean.name = 'IC50_Cleaned'
    Inhibition_clean.name = 'INH_Cleaned'

    # Output csv files
    IC50_clean.to_csv(os.path.join(path, 'IC50_Cleaned.csv'), index=False)
    Inhibition_clean.to_csv(os.path.join(path, 'INH_Cleaned.csv'), index=False)

    # Output smi files
    tempdf = IC50_clean[['SMILES_desalt', 'molecule_chembl_id']]
    tempdf.to_csv(os.path.join(path, 'IC50_Cleaned.smi'), sep='\t', header=False, index=False)
    tempdf = Inhibition_clean[['SMILES_desalt', 'molecule_chembl_id']]
    tempdf.to_csv(os.path.join(path, 'INH_Cleaned.smi'), sep='\t', header=False, index=False)

    # Output Active vs. Inactive Barchart & Donut
    graphJSON3 = plot_status(IC50_clean, Inhibition_clean, path)
    allgraph.append(graphJSON3)


    with open(os.path.join(path, "_graphpreprocess.txt"), "w") as f:
        for s in allgraph:
            f.write(str(s) + "\n")

    response = 'ok'
    return response

def getgraphbioact(filename, uniID):

    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    graphList = []
    with open(os.path.join(path, "_bioact.txt"), "r") as f:
        for line in f:
            graphList.append(line)

    return graphList


def getgraphpreprocess(filename, uniID):

    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    graphList = []
    with open(os.path.join(path, "_graphpreprocess.txt"), "r") as f:
        for line in f:
            graphList.append(line)

    return graphList


def RO5(filename, uniID):
    graphro5 = []

    path = os.path.join(server.config['UPLOAD_FOLDER'], filename+"_"+uniID)
    df = pd.read_csv(os.path.join(path, "IC50_Cleaned.csv"), header=0)

    ID = df.molecule_chembl_id
    pIC50 = df.pIC50
    smiles = df.SMILES_desalt
    #len(smiles)

    from rdkit.Chem import Descriptors
    from rdkit.Chem import MolFromSmiles

    mols = []

    for i in df.SMILES_desalt:
        mol = MolFromSmiles(i)
        mols.append(mol)
    MW = [Descriptors.MolWt(n) for n in mols]
    LogP = [Descriptors.MolLogP(o) for o in mols]
    nHAcc = [Descriptors.NumHAcceptors(p) for p in mols]
    nHDon = [Descriptors.NumHDonors(q) for q in mols]
    print(df['molecule_structures'])
    stdInChiKey = [ast.literal_eval(k)['standard_inchi_key'] if str(k) != "nan" else "-" for k in df['molecule_structures']]
    formula = [ast.literal_eval(k)['full_molformula'] if str(k) != "nan" else "-" for k in df['molecule_properties']]
    smiles = df['canonical_smiles']
    preferredCompoundName = []
    for i in df['molecule_pref_name']:
        if i is not np.NaN:
            preferredCompoundName.append(i)
        else:
            preferredCompoundName.append('-')

    data = pd.DataFrame(
        {'molecule_chembl_id': ID,
         'stdInChiKey': stdInChiKey,
         'formula': formula,
         'STATUS': df.STATUS,
         'pIC50': pIC50,
         'MW': MW,
         'LogP': LogP,
         'nHAcc': nHAcc,
         'nHDon': nHDon,
         'smiles': smiles,
         'preferredCompoundName': preferredCompoundName
         })
    data = data[['molecule_chembl_id', 'stdInChiKey', 'formula', 'STATUS', 'pIC50', 'MW', 'LogP', 'nHAcc', 'nHDon', 'smiles', 'preferredCompoundName']]


    data.to_csv(os.path.join(path, 'IC50_RO5.csv'), sep=',', index=False)

    fig = make_subplots(rows=1, cols=4,
                        subplot_titles=('Molecular Weight',
                                        'Octanol-water Partition<br>Coefficient (logP)',
                                        'Hydrogen Bond Acceptors', 'Hydrogen Bond Donors'),
                        shared_yaxes=True)
    fig.add_trace(
        go.Scatter(x=data[data.STATUS == 'active'].MW, y=data[data.STATUS == 'active'].pIC50, mode='markers',
                   marker_color='#ff3030', name='active'),
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(x=data[data.STATUS == 'inactive'].MW, y=data[data.STATUS == 'inactive'].pIC50, mode='markers',
                   marker_color='#00bfff', name='inactive'),
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(x=data[data.STATUS == 'active'].LogP, y=data[data.STATUS == 'active'].pIC50, mode='markers',
                   marker_color='#ff3030', name='active'),
        row=1, col=2
    )
    fig.add_trace(
        go.Scatter(x=data[data.STATUS == 'inactive'].LogP, y=data[data.STATUS == 'inactive'].pIC50, mode='markers',
                   marker_color='#00bfff', name='inactive'),
        row=1, col=2
    )
    fig.add_trace(
        go.Scatter(x=data[data.STATUS == 'active'].nHAcc, y=data[data.STATUS == 'active'].pIC50, mode='markers',
                   marker_color='#ff3030', name='active'),
        row=1, col=3
    )
    fig.add_trace(
        go.Scatter(x=data[data.STATUS == 'inactive'].nHAcc, y=data[data.STATUS == 'inactive'].pIC50, mode='markers',
                   marker_color='#00bfff', name='inactive'),
        row=1, col=3
    )

    fig.add_trace(
        go.Scatter(x=data[data.STATUS == 'active'].nHDon, y=data[data.STATUS == 'active'].pIC50, mode='markers',
                   marker_color='#ff3030', name='active'),
        row=1, col=4
    )
    fig.add_trace(
        go.Scatter(x=data[data.STATUS == 'inactive'].nHDon, y=data[data.STATUS == 'inactive'].pIC50, mode='markers',
                   marker_color='#00bfff', name='inactive'),
        row=1, col=4
    )
    fig['layout']['xaxis1'].update(title='MW (daltons)')
    fig['layout']['xaxis2'].update(title='logP')
    fig['layout']['xaxis3'].update(title='nHAcc')
    fig['layout']['xaxis4'].update(title='nHDon')
    fig.update_layout(height=500, showlegend=False, yaxis_title='pIC50',margin=dict(l=10, r=10, t=60, b=10))

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    with open(os.path.join(path, "_scatterro5.txt"), 'w') as outfile:
        json.dump(graphJSON, outfile)

    (smw, emw) = defineoutlier_modifiedzscore(data.MW)
    smw = 0
    if emw < 500: emw = 500
    (slp, elp) = defineoutlier_modifiedzscore(data.LogP)
    slp = 0
    if elp <= 5: elp = 5
    (sha, eha) = defineoutlier_modifiedzscore(data.nHAcc)
    sha = 0
    if eha <= 10: eha = 10
    (shd, ehd) = defineoutlier_modifiedzscore(data.nHDon)
    shd = 0
    if ehd <= 5: ehd = 5

    ro5value = ['MW', 'LogP', 'nHAcc', 'nHDon']
    cutoffvalue = ['500', '5', '10',  '5']
    cutoffvalue_start = [smw, slp, sha, shd]
    cutoffvalue_end = [emw, elp, eha, ehd]
    typecom = ['active', 'inactive']

    for ro in ro5value:

        df = data
        # df.set_index('molecule_chembl_id', inplace=True)

        colors = {'active': '#ff3030',
                  'inactive': '#00bfff'}

        fig = go.Figure()

        for i in range(len(typecom)):
            fig.add_trace(go.Histogram(
                x=df[df["STATUS"]==typecom[i]][ro],
                marker_color=colors[typecom[i]],
                opacity=0.6,
                name=typecom[i])
            )

        cutoff = cutoffvalue[ro5value.index(ro)]

        fig.add_shape(
            go.layout.Shape(
                type='line',
                xref='x', yref='paper',
                y0=0, y1=1,
                x0=cutoff, x1=cutoff,
                line=dict(
                    color="gray",
                    width=5,
                    dash="dash"
                )
            )
        )
        mycutoff_start = cutoffvalue_start[ro5value.index(ro)]
        fig.add_shape(
            go.layout.Shape(
                type='line',
                xref='x', yref='paper',
                y0=0, y1=1,
                x0=mycutoff_start, x1=mycutoff_start,
                line=dict(
                    color="black",
                    width=1,
                    dash="dot"
                )
            )
        )
        mycutoff_end = cutoffvalue_end[ro5value.index(ro)]
        fig.add_shape(
            go.layout.Shape(
                type='line',
                xref='x', yref='paper',
                y0=0, y1=1,
                x0=mycutoff_end, x1=mycutoff_end,
                line=dict(
                    color="black",
                    width=1,
                    dash="dot"
                )
            )
        )
        fig.update_layout(
            barmode='overlay',
            margin=dict(l=10, r=10, t=30, b=10),
            height=300
        )
        fig.update_traces(opacity=0.5)
        fig.update_xaxes(title_text=ro)
        fig.update_yaxes(title_text='Frequency')

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        graphro5.append(graphJSON)

    with open(os.path.join(path, "_ro5cutoff.txt"), 'w') as f:
        for s in graphro5:
            f.write(str(s) + "\n")


    # Grab DataFrame rows where column has certain values
    data = data[data.MW.between(smw, emw, inclusive=True) &
                data.LogP.between(slp, elp, inclusive=True) &
                data.nHAcc.between(sha, eha, inclusive=True) &
                data.nHDon.between(shd, ehd, inclusive=True)]  # The inclusive (True: <=, False: <)

    data.to_csv(os.path.join(path, 'IC50_RO5_filter.csv'), sep=',', index=False)


    graphList = []
    with open(os.path.join(path, "_ro5cutoff.txt"), "r") as f:
        for line in f:
            graphList.append(line)

    return graphList

def getscatterro5(filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    graphJSON = ""
    with open(os.path.join(path, "_scatterro5.txt")) as json_file:
        graphJSON = json.load(json_file)

    return graphJSON

def defineoutlier_modifiedzscore(ys):
    threshold = 3.5
    median_y = np.median(ys)
    median_absolute_deviation_y = np.median([np.abs(y - median_y) for y in ys])
    start = ((-threshold*median_absolute_deviation_y)/0.6475)+median_y
    end = ((threshold*median_absolute_deviation_y)/0.6475)+median_y
    return (start, end)

def padel(filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename+"_"+uniID)
    padeldescriptor(mol_dir=(os.path.join(path, 'IC50_Cleaned.smi')), d_file=(os.path.join(path, 'PubchemFingerprinter_IC50.csv')), fingerprints=True, retainorder=True, maxruntime=10000)
    padeldescriptor(mol_dir=(os.path.join(path, 'INH_Cleaned.smi')), d_file=(os.path.join(path, 'PubchemFingerprinter_INH.csv')), fingerprints=True, retainorder=True, maxruntime=10000)


def appendStatus(path):
    Status_data = pd.read_csv(os.path.join(path, 'IC50_RO5_filter.csv'), header=0)
    Status_data = Status_data[['molecule_chembl_id', 'STATUS']]

    for Fingerprints in glob.glob(os.path.join(path, 'PubchemFingerprinter_*')):
        Fp = pd.read_csv(Fingerprints)
        Fp_name = os.path.basename(Fingerprints)

        # Fp_name = Fp.replace(".csv", "")

        #         print ('')
        #         print ('In %s dataset...' % (Fp_name))
        #         print ("There are %d status and %d fingerprint." % (len(Status_data), len(Fp)))

        Fp = Fp.rename(columns={'Name': 'molecule_chembl_id'})
        raw = Status_data.merge(Fp, on='molecule_chembl_id', how='inner')

        # Transform STATUS
        sta = []
        for i in raw.STATUS:
            if i == 'active':
                sta.append(1)  # active
            elif i == 'inactive':
                sta.append(-1)  # inactive
            elif i == 'intermediate':
                sta.append(0)  # intermediate
        raw = raw.drop('STATUS', 1)
        raw["Status"] = sta

        raw.to_csv(os.path.join(path, Fp_name), sep=',', index=False)

def plot_roc_curve(path, y_test, preds):
    fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
    roc_auc = metrics.auc(fpr, tpr)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=fpr, y=tpr,
                             mode='lines',
                             name='Random Forest (AUC='+str(roc_auc)[:str(roc_auc).find('.')+3]+')',
                             line=dict(color='#ef6c00')))
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1],
                             mode='lines',
                             name='Reference',
                             line=dict(color='darkgrey')))
    fig.update_layout(xaxis_title='False Positive Rate (FPR)',
                      yaxis_title='True Positive Rate (TPR)',
                      margin=dict(l=10, r=10, t=30, b=10),
                      height=400,
                      legend=dict(x=.2, y=1.3))

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    with open(os.path.join(path, "_plotroccurve.txt"), 'w') as outfile:
        json.dump(graphJSON, outfile)

    # plt.title('Receiver Operating Characteristic')
    # plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
    # plt.legend(loc = 'lower right')
    # plt.plot([0, 1], [0, 1],'r--')
    # plt.xlim([-0.01, 1.01])
    # plt.ylim([-0.01, 1.01])
    # plt.ylabel('True Positive Rate')
    # plt.xlabel('False Positive Rate')
    # plt.show()

def getplotroccurve (filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    graphJSON = ""
    with open(os.path.join(path, "_plotroccurve.txt")) as json_file:
        graphJSON = json.load(json_file)

    return graphJSON

def fit_ml_algo(algo, X_train, y_train, X_test, y_test, cv, path):

    confusionplot = []

    # One Pass
    model = algo.fit(X_train, y_train)
    train_pred = model.predict(X_train)
    test_pred = model.predict(X_test)
    if (isinstance(algo, (RandomForestClassifier))):
        probs = model.predict_proba(X_test)[:, 1]
    else:
        probs = "Not Available"
    acc = round(model.score(X_train, y_train) * 100, 2)

    # CV
    CV_pred = model_selection.cross_val_predict(algo,
                                                X_train,
                                                y_train,
                                                cv=cv,
                                                n_jobs=-1)
    acc_cv = round(metrics.accuracy_score(y_train, CV_pred) * 100, 2)

    acc_ext = round(model.score(X_test, y_test) * 100, 2)

    conf_mat_Train = confusion_matrix(y_train, train_pred)
    conf_mat_CV = confusion_matrix(y_train, CV_pred)
    conf_mat_Test = confusion_matrix(y_test, test_pred)

    print('*******************************************************************')
    print("Confusion matrix Training set:")
    print(conf_mat_Train)
    tn, fp, fn, tp = conf_mat_Train.ravel()
    Sensitivity = round(float(tp) / (float(tp) + float(fn)) * 100, 2)
    Specificity = round(float(tn) / (float(tn) + float(fp)) * 100, 2)

    print('Sensitivity Training set : %0.2f' % Sensitivity)
    print('Specificity Training set : %0.2f' % Specificity)


    z = [[fn, tn],
         [tp, fp]]
    x = ['Positive', 'Negative']
    y = ['Negative', 'Positive']
    z_text = [[str(fn)+'<br>FN', str(tn)+'<br>TN'],
              [str(tp)+'<br>TP', str(fp)+'<br>FP']]
    fig = ff.create_annotated_heatmap(z, x=x, y=y, annotation_text=z_text, colorscale=['#f8f9fa', '#4db6ac'])
    fig.layout.update({'height':300, 'margin':{'l':10, 'r':10, 't':50, 'b':80}})
    fig.layout.yaxis.update({'title': 'Actual'})
    fig.layout.xaxis.update({'title': 'Predicted'})
    fig.layout.annotations = [{'x': 0.5, 'y': -0.4, 'text': '<b>Confusion matrix Training</b><br>Sensitivity: '+str(
        Sensitivity)+'<br>Specificity: '+ str(Specificity), 'showarrow': False, 'xref': "paper", 'yref': "paper"},
                              {'x':'Positive', 'y':'Negative', 'text': 'FN<br>'+str(fn), 'showarrow':False},
                              {'x':'Positive', 'y':'Positive', 'text': 'TP<br>'+str(tp), 'showarrow':False},
                              {'x':'Negative', 'y':'Negative', 'text': 'TN<br>'+str(tn), 'showarrow':False},
                              {'x':'Negative', 'y':'Positive', 'text': 'FP<br>'+str(fp), 'showarrow':False}]
    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    confusionplot.append(graphJSON)

    print("Confusion matrix CV:")
    print(conf_mat_CV)
    tn, fp, fn, tp = conf_mat_CV.ravel()
    SensitivityCV = round(float(tp) / (float(tp) + float(fn)) * 100, 2)
    SpecificityCV = round(float(tn) / (float(tn) + float(fp)) * 100, 2)

    print('Sensitivity CV : %0.2f' % SensitivityCV)
    print('Specificity CV : %0.2f' % SpecificityCV)

    z = [[fn, tn],
         [tp, fp]]
    x = ['Positive', 'Negative']
    y = ['Negative', 'Positive']
    z_text = [[str(fn) + '<br>FN', str(tn) + '<br>TN'],
              [str(tp) + '<br>TP', str(fp) + '<br>FP']]
    fig = ff.create_annotated_heatmap(z, x=x, y=y, annotation_text=z_text, colorscale=['#f8f9fa', '#4db6ac'])
    fig.layout.update({'height':300, 'margin':{'l':10, 'r':10, 't':50, 'b':80}})
    fig.layout.yaxis.update({'title': 'Actual'})
    fig.layout.xaxis.update({'title': 'Predicted'})
    fig.layout.annotations = [{'x': 0.5, 'y': -0.4, 'text': '<b>Confusion matrix CV</b><br>Sensitivity: '+str(
        SensitivityCV)+'<br>Specificity: '+ str(SpecificityCV), 'showarrow': False, 'xref': "paper", 'yref': "paper"},
                              {'x':'Positive', 'y':'Negative', 'text': 'FN<br>'+str(fn), 'showarrow':False},
                              {'x':'Positive', 'y':'Positive', 'text': 'TP<br>'+str(tp), 'showarrow':False},
                              {'x':'Negative', 'y':'Negative', 'text': 'TN<br>'+str(tn), 'showarrow':False},
                              {'x':'Negative', 'y':'Positive', 'text': 'FP<br>'+str(fp), 'showarrow':False}]
    graphJSON2 = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    confusionplot.append(graphJSON2)

    print("Confusion matrix Testing set:")
    print(conf_mat_Test)
    tn, fp, fn, tp = conf_mat_Test.ravel()
    SensitivityTr = round(float(tp) / (float(tp) + float(fn)) * 100, 2)
    SpecificityTr = round(float(tn) / (float(tn) + float(fp)) * 100, 2)

    print('Sensitivity Testing set : %0.2f' % SensitivityTr)
    print('Specificity Testing set : %0.2f' % SpecificityTr)

    z = [[fn, tn],
         [tp, fp]]
    x = ['Positive', 'Negative']
    y = ['Negative', 'Positive']
    z_text = [[str(fn) + '<br>FN', str(tn) + '<br>TN'],
              [str(tp) + '<br>TP', str(fp) + '<br>FP']]
    fig = ff.create_annotated_heatmap(z, x=x, y=y, annotation_text=z_text, colorscale=['#f8f9fa', '#4db6ac'])
    fig.layout.update({'height': 300, 'margin': {'l': 10, 'r': 10, 't': 50, 'b': 80}})
    fig.layout.yaxis.update({'title': 'Actual'})
    fig.layout.xaxis.update({'title': 'Predicted'})
    fig.layout.annotations = [{'x': 0.5, 'y': -0.4, 'text': '<b>Confusion matrix Testing</b><br>Sensitivity: ' + str(
        SensitivityTr) + '<br>Specificity: ' + str(SpecificityTr), 'showarrow': False, 'xref': "paper", 'yref': "paper"},
                              {'x': 'Positive', 'y': 'Negative', 'text': 'FN<br>'+str(fn), 'showarrow': False},
                              {'x': 'Positive', 'y': 'Positive', 'text': 'TP<br>'+str(tp), 'showarrow': False},
                              {'x': 'Negative', 'y': 'Negative', 'text': 'TN<br>'+str(tn), 'showarrow': False},
                              {'x': 'Negative', 'y': 'Positive', 'text': 'FP<br>'+str(fp), 'showarrow': False}]
    graphJSON3 = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    confusionplot.append(graphJSON3)

    print('')
    # Matthews correlation coefficient (MCC)
    mcc_tr = round(metrics.matthews_corrcoef(y_train, train_pred), 2)
    mcc_cv = round(metrics.matthews_corrcoef(y_train, CV_pred), 2)
    mcc_ext = round(metrics.matthews_corrcoef(y_test, test_pred), 2)

    CV_scores = cross_val_score(algo,
                                X_train,
                                y_train,
                                cv=cv,
                                n_jobs=-1)
    print("Accuracy CV: %0.2f (+/- %0.2f)" % (CV_scores.mean(), CV_scores.std() * 2))

    with open(os.path.join(path, "_confusion.txt"), 'w') as f:
        for s in confusionplot:
            f.write(str(s) + "\n")

    return (CV_pred, test_pred,
            acc, acc_cv, acc_ext, probs,
            mcc_tr, mcc_cv, mcc_ext,
            Sensitivity, Specificity,
            SensitivityCV, SpecificityCV,
            SensitivityTr, SpecificityTr)

def getconfusion (filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    graphList = []
    with open(os.path.join(path, "_confusion.txt"), "r") as f:
        for line in f:
            graphList.append(line)

    return graphList


def Remove_useless_descriptor(data, threshold):
    data2 = (data.drop('Status', axis=1))

    des1 = len(data2.columns)

    h = data2.columns.tolist()
    df = data2.values.astype(np.float)
    df = np.array(df)

    STDEV = np.std(df, axis=0)
    idx = [idx for idx, val in enumerate(STDEV) if val > threshold]
    df2 = df[:, idx]
    hx = np.array(h)[idx]

    df = pd.DataFrame(df2, columns=[hx])

    des2 = len(df.columns)

    df['Status'] = data['Status']

    #     print ('from Remove useless descriptor...')
    #     print ("The initial set of %d descriptors has been reduced to %d descriptors." % (des1, des2))

    return df, des1, des2


def correlation(data, threshold):
    data.columns = data.columns.get_level_values(0)  # Flatten the DF to 2 dimensions
    temp = data.copy()
    des3 = len((data.drop('Status', axis=1).columns))
    corr = stats.pearsonr
    col_corr = set()  # Set of all the names of deleted columns
    corr_matrix = data.corr()
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if corr_matrix.iloc[i, j] >= threshold:
                colname = corr_matrix.columns[i]  # getting the name of column
                col_corr.add(colname)
                if colname in data.columns:
                    del data[colname]  # deleting the column from the dataset
    data['Status'] = temp['Status'].values
    des4 = len((data.drop('Status', axis=1).columns))

    #     print ('from Remove correlation...')
    #     print ("The initial set of %d descriptors has been reduced to %d descriptors." % (des3, des4))

    return data, des3, des4


def STATUS_pca(df):
    STATUS = []

    for i in df.Status:
        if i == 1:
            STATUS.append("active")  # active

        elif i == -1:
            STATUS.append("inactive")  # inactive

    df_pca = df.drop('Status', 1)
    df_pca['STATUS'] = STATUS

    return df_pca


def PCA_variance(data, Fp_name, path):
    # Calculating PCA for both datasets, and graphing the Variance for each feature, per dataset
    from sklearn import preprocessing
    from sklearn.decomposition import PCA

    try:
        std_scale = preprocessing.StandardScaler().fit(data.drop('Status', axis=1))
        X = std_scale.transform(data.drop('Status', axis=1))
        pca1 = PCA(n_components=len(data.columns) - 1)
        fit1 = pca1.fit(X)

        # Graphing the variance per feature
        plt.style.use('seaborn-whitegrid')
        plt.figure(figsize=(25, 7))

        plt.xlabel('PCA Feature')
        plt.ylabel('Variance')
        plt.title('PCA for Discretised Dataset')
        plt.bar(range(0, fit1.explained_variance_ratio_.size), fit1.explained_variance_ratio_);
        plt.savefig(os.path.join(path, Fp_name+'_Variance.pdf'), dpi=300)
    except:
        # removing directory
        shutil.rmtree(path)
        print(path + " doesn't contain enough samples.")
        print("Removing the dataset")
        return -1
    return 1


def PCA_plot(data, Fp_name, path):
    # PCA's components graphed in 2D and 3D
    # Apply Scaling
    from sklearn.decomposition import PCA
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn.preprocessing import StandardScaler
    from matplotlib import pyplot as plt
    data_pca = STATUS_pca(data)

    try:
        # Apply Scaling
        X = data_pca.drop('STATUS', axis=1).as_matrix().astype(np.float)
        y = data_pca['STATUS'].values

        # Formatting
        target_names = ['inactive', 'active']
        colors = ['red', 'blue']
        lw = 2
        alpha = 0.3

        # 2 Components PCA
        plt.style.use('seaborn-whitegrid')
        plt.figure(2, figsize=(20, 8))

        plt.subplot(1, 2, 1)
        pca = PCA(n_components=2)
        X_std = StandardScaler().fit_transform(X)
        #         print(X_std)
        X_r = pca.fit_transform(X_std)
        # X_r = pca.fit(X).transform(X)

        for color, i, target_name in zip(colors, ['inactive', 'active'], target_names):
            plt.scatter(X_r[y == i, 0], X_r[y == i, 1],
                        color=color,
                        alpha=alpha,
                        lw=lw,
                        label=target_name)
        plt.legend(loc='best', shadow=False, scatterpoints=1)
        plt.title('First two PCA directions')

        # 3 Components PCA
        ax = plt.subplot(1, 2, 2, projection='3d')

        pca = PCA(n_components=3)
        X_reduced = pca.fit_transform(X_std)
        for color, i, target_name in zip(colors, ['inactive', 'active'], target_names):
            ax.scatter(X_reduced[y == i, 0], X_reduced[y == i, 1], X_reduced[y == i, 2],
                       color=color,
                       alpha=alpha,
                       lw=lw,
                       label=target_name)
        plt.legend(loc='best', shadow=False, scatterpoints=1)
        ax.set_title("First three PCA directions")
        ax.set_xlabel("1st eigenvector")
        ax.set_ylabel("2nd eigenvector")
        ax.set_zlabel("3rd eigenvector")

        # rotate the axes
        ax.view_init(30, 10)
        plt.savefig(os.path.join(path, Fp_name+'_PCA.pdf'), dpi=300)
    except:
        # removing directory
        shutil.rmtree(path)
        print(path + " doesn't contain enough samples.")
        print("Removing the dataset")
        return -1
    return 1


# Calculating RFE for non-discretised dataset, and graphing the Importance for each feature, per dataset
def Feature_Elimination(data, Fp_name, path):
    start_time = time.time()

    selector1 = RFECV(LogisticRegression(max_iter=1000), step=1, cv=5, n_jobs=-1)
    selector1 = selector1.fit(data.drop('Status', axis=1).values, data['Status'].values)

    print("Feature Ranking For Non-Discretised: %s" % selector1.ranking_)
    print("Optimal number of features : %d" % selector1.n_features_)
    # Plot number of features VS. cross-validation scores
    # plt.style.use('seaborn-whitegrid')
    # plt.figure(figsize=(20, 5))
    # plt.xlabel("Number of features selected - Non-Discretised")
    # plt.ylabel("Cross validation score (nb of correct classifications)")
    # plt.plot(range(1, len(selector1.grid_scores_) + 1), selector1.grid_scores_)
    # plt.savefig(os.path.join(path, 'Feature Elimination.pdf'), dpi=300)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(range(1, len(selector1.grid_scores_) + 1)), y=selector1.grid_scores_,
                             mode='lines',
                             name='Cross validation score',
                             line=dict(color='#4db6ac')))
    fig.update_layout(xaxis_title='Number of features selected - Non-Discretised',
                      yaxis_title='Cross validation score (number of corrects)',
                      margin=dict(l=10, r=10, t=30, b=10),
                      height=350)

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    with open(os.path.join(path, "_elimination.txt"), 'w') as outfile:
        json.dump(graphJSON, outfile)

    # Feature space could be subsetted like so:
    data_ = data.drop('Status', axis=1)
    data2 = data[data.columns[np.insert(selector1.support_, -1, True)]]
    data2['Status'] = data['Status'].values
    data2.to_csv(os.path.join(path, Fp_name+' Feature Elimination.csv'), sep=',', index=False)

    print("Feature optimization took %.2f seconds"
          % ((time.time() - start_time)))

    return data2, selector1.n_features_, selector1.support_

def getelimination(filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    graphJSON = ""
    with open(os.path.join(path, "_elimination.txt")) as json_file:
        graphJSON = json.load(json_file)

    return graphJSON

def RFC(X_train, y_train, X_test, y_test, path, Fp_name, data2):
    import time, datetime
    from sklearn.ensemble import RandomForestClassifier
    from sklearn import model_selection, metrics
    import pickle

    start_time = time.time()
    rfc = RandomForestClassifier(bootstrap=False, class_weight=None,
                                 criterion='entropy', max_depth=None, max_features='sqrt',
                                 max_leaf_nodes=None, min_impurity_decrease=1e-07,
                                 min_samples_leaf=1, min_samples_split=2,
                                 min_weight_fraction_leaf=0.0, n_estimators=44, n_jobs=1,
                                 oob_score=False, random_state=4, verbose=False,
                                 warm_start=False)

    (train_pred_PC,
     test_pred_PC,
     acc_PC,
     acc_cv_PC,
     acc_ext_PC,
     probs_PC,
     mcc_PC,
     mcc_cv_PC,
     mcc_ext_PC,
     Sen_PC,
     Spe_PC,
     SenCV_PC,
     SpeCV_PC,
     SenTr_PC,
     SpeTr_PC) = fit_ml_algo(rfc, X_train, y_train, X_test, y_test, 10, path)

    # save the model to disk
    # uniID = path.replace('dataset\\', '')
    # uniID = uniID.replace('\\', '')
    filenameX = Fp_name + '_RF_model.sav'
    # Making Directory
    # if not os.path.exists(os.path.join(path, 'model')):
    #     os.makedirs(os.path.join(path, 'model'))
    # newpath = os.path.join(path, 'model')
    # pickle.dump(rfc, open(os.path.join(newpath, '_' +filenameX), 'wb'))
    pickle.dump(rfc, open(os.path.join(path, filenameX), 'wb'))



    PC_time = (time.time() - start_time)

    print('')
    print("Running Time: %s" % datetime.timedelta(seconds=PC_time))
    print("Accuracy: %s" % acc_PC)

    print("Accuracy CV 10-Fold: %s" % acc_cv_PC)

    print("MCC train: %s" % mcc_PC)
    print("MCC CV: %s" % mcc_cv_PC)
    print("MCC Test: %s" % mcc_ext_PC)

    print("Accuracy External: %s" % acc_ext_PC)

    classifyreport = []

    print('Classification report: Training set')
    print(metrics.classification_report(y_train, train_pred_PC))

    classifyreport.append('<b class="stattitle">Classification-Report: Training-set</b><br>')
    metrics_train = str(metrics.classification_report(y_train, train_pred_PC))
    mtrain = metrics_train.split()
    table_train = '<table class="tablet">' \
                  ' <tr>' \
                  '     <th></th><th>'+mtrain[0]+'</th><th>'+mtrain[1]+'</th><th>'+mtrain[2]+'</th><th>'+mtrain[3]+'</th>' \
                  ' </tr>' \
                  ' <tr>' \
                  '     <td>'+mtrain[4]+'</td><td>'+mtrain[5]+'</td><td>'+mtrain[6]+'</td><td>'+mtrain[7]+'</td><td>'+mtrain[8]+'</td>' \
                  ' </tr>' \
                  ' <tr>' \
                  '     <td>'+mtrain[9]+'</td><td>'+mtrain[10]+'</td><td>'+mtrain[11]+'</td><td>'+mtrain[12]+'</td><td>'+mtrain[13]+'</td>' \
                  ' </tr>'\
                  ' <tr>' \
                  '     <td>'+mtrain[14]+'</td><td></td><td></td><td>'+mtrain[15]+'</td><td>'+mtrain[16]+'</td>' \
                  ' </tr>' \
                  ' <tr>' \
                  '     <td>'+mtrain[17]+' '+mtrain[18]+'</td><td>'+mtrain[19]+'</td><td>'+mtrain[20]+'</td><td>'+mtrain[21]+'</td><td>'+mtrain[22]+'</td>' \
                  ' </tr>' \
                  ' <tr>' \
                  '     <td>'+mtrain[23]+' '+mtrain[24]+'</td><td>'+mtrain[25]+'</td><td>'+mtrain[26]+'</td><td>'+mtrain[27]+'</td><td>'+mtrain[28]+'</td>' \
                  ' </tr>' \
                  '</table>'
    classifyreport.append(table_train)

    print('Classification report: Test set')
    print(metrics.classification_report(y_test, test_pred_PC))

    classifyreport.append('<b class="stattitle">Classification-Report: Test-set</b><br>')
    metrics_test = str(metrics.classification_report(y_test, test_pred_PC))
    mtest = metrics_test.split()
    table_test = '<table class="tablet">' \
                  ' <tr>' \
                  '     <th></th><th>'+mtest[0]+'</th><th>'+mtest[1]+'</th><th>'+mtest[2]+'</th><th>'+mtest[3]+'</th>' \
                  ' </tr>' \
                  ' <tr>' \
                  '     <td>'+mtest[4]+'</td><td>'+mtest[5]+'</td><td>'+mtest[6]+'</td><td>'+mtest[7]+'</td><td>'+mtest[8]+'</td>' \
                  ' </tr>' \
                  ' <tr>' \
                  '     <td>'+mtest[9]+'</td><td>'+mtest[10]+'</td><td>'+mtest[11]+'</td><td>'+mtest[12]+'</td><td>'+mtest[13]+'</td>' \
                  ' </tr>'\
                  ' <tr>' \
                  '     <td>'+mtest[14]+'</td><td></td><td></td><td>'+mtest[15]+'</td><td>'+mtest[16]+'</td>' \
                  ' </tr>' \
                  ' <tr>' \
                  '     <td>'+mtest[17]+' '+mtest[18]+'</td><td>'+mtest[19]+'</td><td>'+mtest[20]+'</td><td>'+mtest[21]+'</td><td>'+mtest[22]+'</td>' \
                  ' </tr>' \
                  ' <tr>' \
                  '     <td>'+mtest[23]+' '+mtest[24]+'</td><td>'+mtest[25]+'</td><td>'+mtest[26]+'</td><td>'+mtest[27]+'</td><td>'+mtest[28]+'</td>' \
                  ' </tr>' \
                  '</table>'
    classifyreport.append(table_test)

    with open(os.path.join(path, "_statvalue.txt"), "w") as f:
        for s in classifyreport:
            f.write(str(s) + "\n")


    plot_roc_curve(path, y_test, probs_PC)
    # metrics.plot_roc_curve(rfc, X_test, y_test)

    #     Feature_All_PC     = des1
    #     Feature_STDEV_PC   = des2
    #     Faature_Corr_PC    = des4
    #     Feature_Seclect_PC = 59

    # Using Random Forest to gain an insight on Feature Importance
    plt.style.use('seaborn-whitegrid')
    # feature_importance = clf.feature_importances_
    importances = rfc.feature_importances_
    importance = pd.DataFrame(importances, index=data2.drop('Status', axis=1).columns, columns=["Importance"])

    importance = importance.sort_values(by='Importance', ascending=True).tail(10)
    # importance.plot(kind='barh', figsize=(12, len(importance) / 2))  #
    # plt.tick_params(axis='both', which='major', labelsize=14)
    # plt.xlabel('Feature importance', fontsize=16, fontweight='bold')
    # plt.gca().legend_.remove()
    # # plt.gca().tick_params(axis='both', which='major', pad=15)
    # plt.savefig(os.path.join(path, Fp_name+' Importance_RF.pdf'), dpi=300)

    fig = go.Figure(go.Bar(
        x=list(importance['Importance']),
        y=list(importance.index),
        orientation='h',
        text=list(importance.index),
        textposition="outside",
        name="Feature importance",
        marker_color='#ef6c00')
    )
    fig.update_layout(
        autosize=True,
        margin=dict(l=10, r=10, t=30, b=10),
        showlegend=False,
        height=350
    )
    fig.update_xaxes(title_text='Feature Importance of Random Forest Classifier')
    fig.update_yaxes(showticklabels=False)
    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    with open(os.path.join(path, "_importance.txt"), 'w') as outfile:
        json.dump(graphJSON, outfile)


def getstatvalue(filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    Row_list = ''
    with open(os.path.join(path, "_statvalue.txt"), "r") as f:
        for line in f:
            Row_list += line.replace('\n', '<br>')

    return Row_list

def getimportance(filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    graphJSON = ''
    with open(os.path.join(path, "_importance.txt")) as json_file:
        graphJSON = json.load(json_file)

    return graphJSON


# split into training and test data using 80/20 ratio
sss = StratifiedShuffleSplit(n_splits=10, test_size=0.2, random_state=0)


def Balance(data2, data, Fp_name, path):
    X = (data2.drop('Status', axis=1)).values
    y = data2['Status'].values

    sss.get_n_splits(X, y)

    for train_index, test_index in sss.split(X, y):
        #         print("TRAIN:", train_index, "TEST:", test_index)
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

    # Instanciate a PCA object for the sake of easy visualisation
    pca = PCA(n_components=2)
    # Fit and transform x to visualise inside a 2D feature space
    X_vis = pca.fit_transform(X_train)

    # Apply SMOTE + Tomek links
    sm = SMOTEENN(random_state=0)
    X_resampled, y_resampled = sm.fit_sample(X_train, y_train)
    if len(X_resampled) == 0 or len(y_resampled) == 0:
        # removing directory
        shutil.rmtree(path)
        print(path + " doesn't contain enough samples.")
        print("Removing the dataset")
        return -1, -1, -1, -1, -1
    X_res_vis = pca.transform(X_resampled)

    print("The dataset present %d compounds and train %d compounds has been resampled to %d compounds." % (
    len(y), len(y_train), len(y_resampled)))

    return X_resampled, y_resampled, X_test, y_test, 1


def modeling(filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename+"_"+uniID)

    appendStatus(path)
    final_time = time.time()
    Fingerprints = (os.path.join(path, "PubchemFingerprinter_IC50.csv"))

    data = pd.read_csv(Fingerprints)

    Fp_name = os.path.basename(Fingerprints)
    Fp_name = Fp_name.replace(".csv", "")
    Fp_name = Fp_name.replace("PubchemFingerprinter_", "")

    data = data.iloc[:, 1:]
    data = data.drop(data.index[data.Status == 0])
    data.index = range(len(data))

    data, des1, des2 = Remove_useless_descriptor(data, 0.05)  # Remove correlation cut off 95%
    data, des3, des4 = correlation(data, 0.7)  # RemFp_name, ove correlation cut off 0.7

    data.columns = data.columns.get_level_values(0)  # Flatten the DF to 2 dimensions

    # Using Random Forest to gain an insight on Feature Importance
    clf = RandomForestClassifier()
    clf.fit(data.drop('Status', axis=1), data['Status'])

    plt.style.use('seaborn-whitegrid')
    # feature_importance = clf.feature_importances_
    importances = clf.feature_importances_
    importance = pd.DataFrame(importances, index=data.drop('Status', axis=1).columns, columns=["Importance"])

    importance = importance.sort_values(by='Importance', ascending=True).tail(20)
    # importance.plot(kind='barh', figsize=(10, len(importance) / 2))  #
    # plt.savefig(os.path.join(path, Fp_name + '_Feature_Importance.pdf'), dpi=300)

    data2, Feature, support = Feature_Elimination(data, Fp_name, path)

    #     print(uniID)
    # cp1 = PCA_variance(data2, Fp_name, path)
    # cp2 = PCA_plot(data2, Fp_name, path)
    # if cp1 == -1 or cp2 == -1:
    #     print("Skipping to the next dataset")
    #     # return

    # Balance and split sample
    X_train, y_train, X_test, y_test, cp = Balance(data2, data, Fp_name, path)
    if cp == -1:
        print("Skipping to the next dataset")
        # return

    # RandomForest
    RFC(X_train, y_train, X_test, y_test, path, Fp_name, data2)

    print("--------------------------------------")

def create_datatablecompound(filename, uniID):

    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)

    df = pd.read_csv(os.path.join(path, "IC50_RO5_filter.csv"))
    html = """<table id="tablecompound_div" class="table display row-border">
        <thead>
            <tr>"""
    i = 0
    for header in df.columns.values:
        html += "<th>" + header + "</th>"
        if i == 0:
            html += "<th>2D</th>"
            i = 1
    html += """</tr>
        </thead>
    </table>"""

    return html


def create_datacompound(filename, uniID):
    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    df = pd.read_csv(os.path.join(path, "IC50_RO5_filter.csv"))

    Row_list = []

    # Iterate over each row
    for index, rows in df.iterrows():
        # append the list to the final list
        templist = list(rows.values)
        # print(templist[9], end=' ')
        from rdkit.Chem import AllChem, Draw
        m = Chem.MolFromSmiles(str(templist[9]))
        fig = Draw.MolToImage(m)
        import base64
        from io import BytesIO
        buff = BytesIO()
        fig.save(buff, format="PNG")
        img_str = base64.b64encode(buff.getvalue())
        templist.insert(1, "<img style='width: 100px; height: 100px' src='data:image/png;base64,"+str(img_str)[2:-1]+"'>")
        Row_list.append(templist)

    with open(os.path.join(path, "_datacompound.txt"), "w") as f:
        for s in Row_list:
            f.write(str(s) + "\n")

    Row_list = []
    with open(os.path.join(path, "_datacompound.txt"), "r") as f:
        for line in f:
            Row_list.append(ast.literal_eval(line))


    return Row_list


def getcompoundinfo(inchikey):
    mc = get_client('chem')
    k = inchikey
    m = mc.query('pubchem.inchi_key:' + k, as_dataframe=True)
    cid = '-'
    if 'pubchem.cid' in m.columns:
        cid = m['pubchem.cid'].values[0]
    iupac = '-'
    if 'pubchem.iupac.traditional' in m.columns:
        iupac = m['pubchem.iupac.traditional'].values[0]
    print(cid, iupac)

    return [str(cid), str(iupac)]


def decodeFP(str):
   fp_hex = base64.b64decode(str).hex()
   return '{0:020b}'.format(int(fp_hex[8:], 16))[:-7].zfill(881)


def predict(filename, uniID, filecompounds):

    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)

    pathfilecompounds = os.path.join(path, 'compoundstest')

    # if os.path.exists(os.path.join(pathfilecompounds, "_"+filecompounds+"_Predicted.csv")):
    #     compoundstest_result = pd.read_csv(os.path.join(pathfilecompounds, "_" + filecompounds + "_Predicted.csv"))
    #     graph = predictvisualize(pathfilecompounds, compoundstest_result)
    #     return graph

    com = pd.read_csv(os.path.join(pathfilecompounds, filecompounds))
    print(com)

    duplicateRowsDF = com[com.duplicated()]
    com = com.drop_duplicates(subset="InChIKey", keep='first', inplace=False)
    com['InChIKey'].replace('', np.nan, inplace=True)
    com.dropna(subset=['InChIKey'], inplace=True)
    if not duplicateRowsDF.empty:
        print("Duplicate Keys: " + duplicateRowsDF)

    keys = ''
    count = 0
    df = pd.DataFrame()
    for i in tqdm(com['InChIKey']):
        keys = keys + str(i) + ","
        count = count + 1
        if (count % 100) == 0 or count == len(com['InChIKey']):
            keys = keys[:-1]
            # /inchikey,CanonicalSMILES,IsomericSMILES,Fingerprint2D request parameters
            resp = requests.get(
                'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/' + keys + '/property/inchikey,CanonicalSMILES,IsomericSMILES,Fingerprint2D/JSON')
            if resp.status_code != 200:
                # This means something went wrong.
                print("Error: " + str((resp.status_code)))
            keys = ''
            res = resp.json()
            for i in res.values():
                for j in i.values():
                    temp = pd.DataFrame(j)
            df = df.append(temp)

    # fp = pd.DataFrame()
    col_names = []
    # fp['InChIKey'] = df['InChIKey']
    for i in range(881):
        col_name = 'PubchemFP' + str(i)
        col_names.append(col_name)
    #     fp[col_name] = 0
    fp = pd.DataFrame(columns=['InChIKey']+col_names)
    fp_list = []
    for i in df['Fingerprint2D'].tolist():
        fp_list.append(decodeFP(i))
    df['hexFP'] = fp_list
    compound_index = 0
    fp_index = 1
    # for i in tqdm(df['hexFP']):
    simple = []
    for i in tqdm(range(len(df['hexFP']))):
        fp_index = 1
        # print(list(df.iloc[i]['hexFP']))
        simple.append([df.iloc[i]['InChIKey']]+list(df.iloc[i]['hexFP']))
        # for j in list(i):
        #     fp.iloc[compound_index, fp_index] = j
        #     fp_index += 1
        # compound_index += 1
    fp = pd.DataFrame(simple, columns=['InChIKey']+col_names)
    fp.to_csv(os.path.join(pathfilecompounds, filecompounds + '_test.csv'), index=False)

    ### Old Fingerprint Calculation ##
    # df = df.rename(columns={'CanonicalSMILES': 'canonical_smiles', 'IsomericSMILES': 'isomeric_smiles'})
    # print(df.columns)
    # df = clean_smiles(df)
    #
    # df2 = df[['SMILES_desalt', 'InChIKey']]
    # df2.to_csv(os.path.join(pathfilecompounds, filecompounds + '_test.smi'), sep='\t', header=False, index=False)
    #
    # padeldescriptor(mol_dir=os.path.join(pathfilecompounds, filecompounds + '_test.smi'), d_file=os.path.join(pathfilecompounds, filecompounds + '_test.csv'), fingerprints=True, threads=6,
    #                 retainorder=True, maxruntime=10000)

    compoundstest = pd.read_csv(os.path.join(pathfilecompounds, filecompounds + '_test.csv'))
    filename_model = 'IC50_RF_model.sav'
    filename_feature = 'IC50 Feature Elimination.csv'

    pathmodel = os.path.join(path, filename_model)
    loaded_model = pickle.load(open(pathmodel, 'rb'))
    fea_elim = pd.read_csv(os.path.join(path, filename_feature))

    features = []
    for col in fea_elim.columns:
        features.append(col)
    print(features)

    compoundstest = compoundstest.set_index('InChIKey')
    compoundstest = compoundstest.filter(features)
    print(compoundstest)

    y_pred = loaded_model.predict(compoundstest)
    y_prob = loaded_model.predict_proba(compoundstest)

    y_maxprob = []
    for i in y_prob.tolist():
        y_maxprob.append(max(i))

    compoundstest_result = pd.DataFrame()
    compoundstest_result['InChIKey'] = list(compoundstest.index)
    compoundstest_result['CanonicalSMILES'] = list(df['CanonicalSMILES'])
    print(compoundstest_result['CanonicalSMILES'])
    compoundstest_result['Status'] = y_pred.tolist()
    compoundstest_result['Confidence Score'] = y_maxprob

    compoundstest_result.to_csv(os.path.join(pathfilecompounds, "_"+filecompounds+"_Predicted.csv"))

    print(compoundstest_result)
    graph = predictvisualize(pathfilecompounds, compoundstest_result)
    return graph

def predictvisualize(pathfilecompounds, compoundstest_result):

    prediction_active = compoundstest_result[compoundstest_result['Status'] == 1]['Status']
    prediction_inactive = compoundstest_result[compoundstest_result['Status'] == -1]['Status']
    prediction_total = len(prediction_active)+len(prediction_inactive)

    top_labels = ['Active', 'Inactive']
    colors = ['#ff3030', '#00bfff']
    x_data = [[len(prediction_active), len(prediction_inactive)]]
    y_data = ['Compound']

    fig = go.Figure()
    for i in range(0, len(x_data[0])):
        for xd, yd in zip(x_data, y_data):
            fig.add_trace(go.Bar(
                x=[xd[i]], y=[yd],
                orientation='h',
                marker=dict(
                    color=colors[i],
                    line=dict(color='rgb(248, 248, 249)', width=1)
                ),
                name=""
            ))

    fig.update_layout(
        xaxis=dict(
            showgrid=False,
            showline=False,
            showticklabels=False,
            zeroline=False,
            domain=[0.05, 1]
        ),
        yaxis=dict(
            showgrid=False,
            showline=False,
            showticklabels=False,
            zeroline=False
        ),
        barmode='stack',
        paper_bgcolor='#fff',
        plot_bgcolor='#fff',
        margin=dict(l=10, r=10, t=50, b=10),
        showlegend=False
    )

    annotations = []

    for yd, xd in zip(y_data, x_data):

        # labeling the first percentage of each bar (x_axis)
        annotations.append(dict(xref='x', yref='y',
                                x=xd[0] / 2, y=yd,
                                text=str("{0:.2f}".format(xd[0]/prediction_total*100)) + '%',
                                font=dict(size=26,
                                          color='rgb(248, 248, 255)'),
                                showarrow=False))
        # labeling the first Likert scale (on the top)
        if yd == y_data[-1]:
            annotations.append(dict(xref='x', yref='paper',
                                    x=xd[0] / 2, y=1.2,
                                    text=top_labels[0],
                                    font=dict(size=26,
                                              color='rgb(67, 67, 67)'),
                                    showarrow=False))
        space = xd[0]
        for i in range(1, len(xd)):
            # labeling the rest of percentages for each bar (x_axis)
            annotations.append(dict(xref='x', yref='y',
                                    x=space + (xd[i] / 2), y=yd,
                                    text=str("{0:.2f}".format(xd[i]/prediction_total*100)) + '%',
                                    font=dict(size=26,
                                              color='rgb(248, 248, 255)'),
                                    showarrow=False))
            # labeling the Likert scale
            if yd == y_data[-1]:
                annotations.append(dict(xref='x', yref='paper',
                                        x=space + (xd[i] / 2), y=1.2,
                                        text=top_labels[i],
                                        font=dict(size=26,
                                                  color='rgb(67, 67, 67)'),
                                        showarrow=False))
            space += xd[i]

    fig.update_layout(annotations=annotations, height=200)
    fig.update_yaxes(showticklabels=False)

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    with open(os.path.join(pathfilecompounds, "_predictionratio.txt"), 'w') as outfile:
        json.dump(graphJSON, outfile)

    with open(os.path.join(pathfilecompounds, "_predictionratio.txt")) as json_file:
        graphJSON = json.load(json_file)

    return graphJSON


def predictiontable(filename, uniID, filenamecompound):

    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    pathfilecompounds = os.path.join(path, 'compoundstest')

    df = pd.read_csv(os.path.join(pathfilecompounds, "_"+filenamecompound+"_Predicted.csv"))
    df = df.sort_values(by=['Status', 'Confidence Score'], ascending=False)
    df = df.iloc[:, 1:]

    html = """<table id="tablepredict_div" class="table display row-border">
        <thead>
            <tr>"""
    i = 0
    for header in df.columns.values:
        html += "<th>" + header + "</th>"
        if i == 0:
            html += "<th>2D</th>"
            i = 1
    html += """</tr>
        </thead>
    </table>"""

    return html

def create_datapredict(filename, uniID, filenamecompound):

    path = os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + uniID)
    pathfilecompounds = os.path.join(path, 'compoundstest')

    df = pd.read_csv(os.path.join(pathfilecompounds, "_"+filenamecompound+"_Predicted.csv"))
    df = df.sort_values(by=['Status', 'Confidence Score'], ascending=False)
    df = df.iloc[:, 1:]
    df['Status'] = df['Status'].replace({1: 'active', -1: 'inactive'})


    Row_list = []

    # Iterate over each row
    for index, rows in df.iterrows():
        # append the list to the final list
        templist = list(rows.values)
        # print(templist[9], end=' ')
        mc = get_client('chem')
        k = templist[1]
        m = Chem.MolFromSmiles(str(k))
        fig = Draw.MolToImage(m)
        import base64
        from io import BytesIO
        buff = BytesIO()
        fig.save(buff, format="PNG")
        img_str = base64.b64encode(buff.getvalue())
        templist.insert(1, "<img style='width: 100px; height: 100px' src='data:image/png;base64," + str(img_str)[
                                                                                                2:-1] + "'>")
        # templist.insert(1, "No image")
        templist[-1] = "{0:.2f}".format(templist[-1])
        Row_list.append(templist)

    with open(os.path.join(pathfilecompounds, "_"+filenamecompound+"_datapredict.txt"), "w") as f:
        for s in Row_list:
            f.write(str(s) + "\n")

    Row_list = []
    with open(os.path.join(pathfilecompounds, "_"+filenamecompound+"_datapredict.txt"), "r") as f:
        for line in f:
            Row_list.append(ast.literal_eval(line))


    return Row_list


def getcompoundinfopredict(inchikey):
    mc = get_client('chem')
    k = inchikey
    m = mc.query('pubchem.inchi_key:' + k, as_dataframe=True)
    cid = '-'
    if 'pubchem.cid' in m.columns:
        cid = m['pubchem.cid'].values[0]
    iupac = '-'
    if 'pubchem.iupac.traditional' in m.columns:
        iupac = m['pubchem.iupac.traditional'].values[0]
    formula = '-'
    if 'pubchem.molecular_formula' in m.columns:
        formula = m['pubchem.molecular_formula'].values[0]
    smiles = '-'
    if 'pubchem.smiles.canonical' in m.columns:
        smiles = m['pubchem.smiles.canonical'].values[0]
    else:
        smiles = m['pubchem.smiles.isomeric'].values[0]
    prefname = '-'
    if 'chembl.pref_name' in m.columns:
        prefname = m['chembl.pref_name'].values[0]
    chemblid = '-'
    if 'chembl.molecule_chembl_id' in m.columns:
        chemblid = m['chembl.molecule_chembl_id'].values[0]
    return [str(cid), str(iupac), str(formula), str(smiles), str(prefname), str(chemblid)]