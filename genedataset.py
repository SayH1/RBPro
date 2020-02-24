import os
import time

from cffi.backend_ctypes import xrange
from flask import Flask
import pandas as pd
import plotly
import plotly.graph_objs as go
import numpy as np
import json
import requests
import math
import scipy.stats as ss
import ast
from IPython.display import JSON
import xml.etree.ElementTree as ET
from tqdm import tqdm
from pandas import json_normalize

# R Integration
# rpy2 version 2.9.4
import rpy2

print(rpy2.__version__)
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
# import rpy2's package module
import rpy2.robjects.packages as rpackages
from rpy2.robjects import r, pandas2ri
# R vector of strings
from rpy2.robjects.vectors import StrVector

# import R's "base" package
base = importr('base')
# import R's "utils" package
utils = importr('utils')
# import R's utility package
utils = rpackages.importr('utils')
# select a mirror for R packages
# select the first mirror in the list
utils.chooseCRANmirror(ind=1)
# utils.install_packages('BiocManager')
BiocManager = importr('BiocManager')
# BiocManager.install("edgeR")
edgeR = importr('edgeR')
# BiocManager.install("Glimma")
Glimma = importr('Glimma')
# utils.install_packages('dplyr')
dplyr = importr('dplyr')
# utils.install_packages("gplots")
gplots = importr('gplots')
r.source(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'rpackage/rfunction.R'))
# utils.install_packages('tibble')
tibble = importr('tibble')

# Application
dir = os.path.abspath(os.getcwd())
dir = (dir+'\\file\\')
dir = dir.replace('\\', '/')
server = Flask(__name__)
server.config['UPLOAD_FOLDER'] = dir

# Global Variables
data = pd.DataFrame()
geneheader = ''
sampletype = []

# pathway_libraries = ['Reactome_2016']
# pathway_libraries = ['KEGG_2019_Human', 'WikiPathways_2019_Human', 'Reactome_2016']
pathway_libraries = ['KEGG_2019_Human', 'WikiPathways_2019_Human', 'Reactome_2016', 'BioPlanet_2019', 'BioCarta_2016',
                     'Panther_2016']
pathwaylibrariesloop = ['KEGG_2019_Human', 'WikiPathways_2019_Human', 'Reactome_2016']


# Function to get header of gene name and list of sample types
def assign_geneheadersampletype(genehead, sample):
    global geneheader
    geneheader = genehead
    global sampletype
    sampletype = sample


# Function to read the data from the given file name
def read_genedataset(filename):
    data = pd.read_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename))
    data[geneheader] = data[geneheader].astype(str)
    return data


# Function to read the header of data from the given file name
def read_headerdataset(filename):
    data = pd.read_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename))
    headerframe = data.columns.values
    return headerframe


def assign_geneheader(genehead):
    global geneheader
    geneheader = genehead


def read_geneheader():
    return geneheader


def get_pathwaydb():
    return pathway_libraries


def check_exist(filename):
    if os.path.isfile(os.path.join(server.config['UPLOAD_FOLDER'], filename)):
        return True
    else:
        return False


def create_datatable(filename, update):
    print("function: create_datatable()")
    start = time.time()

    df = read_genedataset(filename)
    html = """<table id="table_div" class="table display">
        <thead>
            <tr>"""
    i = 0
    for header in df.columns.values:
        if i == 0:
            html += "<th>" + geneheader + "</th>"
            i += 1
        else:
            html += "<th>" + header + "</th>"
    html += """</tr>
        </thead>
    </table>"""

    end = time.time()
    print("create_datatable", end - start)
    return html


def create_datatablevalues(filename, update):
    print("function: create_datatablevalues()")
    start = time.time()

    Row_list = []

    if update == True:
        df = read_genedataset(filename)
        # Iterate over each row
        for index, rows in df.iterrows():
            # append the list to the final list
            Row_list.append(list(rows.values))

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavalues.txt"), "w") as f:
            for s in Row_list:
                f.write(str(s) + "\n")
    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavalues.txt"), "r") as f:
            for line in f:
                Row_list.append(ast.literal_eval(line))

    end = time.time()
    print("create_datatablevalues", end - start)
    return Row_list


def create_plot(filename, update):
    print("function: create_plot()")
    start = time.time()

    graphJSON = ""

    if update == True:
        data = read_genedataset(filename)

        i = 0
        dataplot = []
        for header in list(data.columns.values):
            if i == 0:
                i += 1
            else:
                dataplot.append(go.Bar(
                    name=header,
                    x=data[geneheader],
                    y=data[header]
                ))
        # Change the bar mode
        fig = go.Figure(data=dataplot)
        fig.update_layout(
            autosize=True,
            barmode='group',
            height=500,
            margin=dict(l=10, r=20, t=30, b=10)
        )
        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_barjson.txt"), 'w') as outfile:
            json.dump(graphJSON, outfile)

    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_barjson.txt")) as json_file:
            graphJSON = json.load(json_file)

    end = time.time()
    print("create_plot", end - start)
    return graphJSON


def create_plotlibrary(filename, update):
    print("function: create_plotlibrary()")
    start = time.time()

    graphJSON = ""

    if update == True:
        data = read_genedataset(filename)
        colors = {'Normal': '#4db6ac',
                  'Perturbation': '#ef6c00'}
        color = [colors[i] for i in sampletype]

        hover = []
        for i in range(len(sampletype)):
            hover.append(sampletype[i] + "<br>" + str(list(data.iloc[:, 1:].sum())[i]))

        fig = go.Figure(go.Bar(
            x=list(data.iloc[:, 1:].sum()),
            y=list(data.columns[1:]),
            orientation='h',
            text=list(data.columns[1:]),
            textposition="inside",
            hovertext=hover,
            name="Library Size",
            marker_color=color)
        )
        fig.update_layout(
            autosize=True,
            margin=dict(l=10, r=10, t=30, b=10),
            showlegend=False
        )
        fig.update_xaxes(title_text='Library Size (Total Read Counts)')
        fig.update_yaxes(showticklabels=False)
        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_libraryjson.txt"), 'w') as outfile:
            json.dump(graphJSON, outfile)

    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_libraryjson.txt")) as json_file:
            graphJSON = json.load(json_file)

    end = time.time()
    print("create_plotlibrary", end - start)
    return graphJSON


def create_histogramcutoff(filename, update):
    print("function: create_histogramcutoff()")
    start = time.time()

    graphJSON = ""

    if update == True:
        df = read_genedataset(filename)
        df.set_index(geneheader, inplace=True)
        lsize = [np.sum(df.iloc[:, i]) for i in range(len(df.columns))]
        data = (df / df.sum()) * (10 ** 6)
        data = data.fillna(0)
        data = np.log2(data + (2 / (np.mean(lsize) / (10 ** 6))))

        colors = {'Normal': '#4db6ac',
                  'Perturbation': '#ef6c00'}

        fig = go.Figure()

        for i in range(len(data.columns)):
            fig.add_trace(go.Histogram(
                x=data.iloc[:, i],
                marker_color=colors[sampletype[i]],
                opacity=0.75,
                name=data.columns[i])
            )

        L = np.mean(lsize) / (10 ** 6)
        M = np.median(lsize) / (10 ** 6)
        cutoff = np.log2(10 / M + 2 / L)

        fig.add_shape(
            go.layout.Shape(
                type='line',
                yref='paper',
                y0=0,
                y1=1,
                xref='x',
                x0=cutoff,
                x1=cutoff,
                line=dict(
                    color="gray",
                    width=1,
                    dash="dot"
                )
            )
        )

        fig.update_layout(
            barmode='overlay',
            margin=dict(l=10, r=10, t=30, b=10)
        )
        fig.update_traces(opacity=0.5)
        fig.update_xaxes(title_text='Log-CPM')
        fig.update_yaxes(title_text='Frequency (Genes)')

        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_histogramcutoff.txt"), 'w') as outfile:
            json.dump(graphJSON, outfile)

    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_histogramcutoff.txt")) as json_file:
            graphJSON = json.load(json_file)

    end = time.time()
    print("create_histogramcutoff", end - start)
    return graphJSON


def create_PCA2D(filename, update):
    print("function: create_PCA2D()")
    start = time.time()

    graphJSON = ""

    if update == True:
        data = read_genedataset(filename)
        data = data.transpose()
        X = data.iloc[1:, :].values
        y = np.array(sampletype)

        from sklearn.preprocessing import StandardScaler
        X_std = StandardScaler().fit_transform(X)

        from sklearn.decomposition import PCA as sklearnPCA
        sklearn_pca = sklearnPCA(n_components=2)
        Y_sklearn = sklearn_pca.fit_transform(X_std)

        dataplot = []
        colors = {'Normal': '#4db6ac',
                  'Perturbation': '#ef6c00'}

        for name, col in zip(('Normal', 'Perturbation'), colors.values()):
            trace = dict(
                type='scatter',
                x=Y_sklearn[y == name, 0],
                y=Y_sklearn[y == name, 1],
                mode='markers',
                name=name,
                marker=dict(
                    color=col,
                    size=12,
                    line=dict(
                        color='rgba(217, 217, 217, 0.14)',
                        width=0.5),
                    opacity=0.8)
            )
            dataplot.append(trace)

        layout = dict(scene=dict(
            xaxis=dict(title='PC1', showline=False),
            yaxis=dict(title='PC2', showline=False)),
            autosize=True,
            height=500,
            margin=dict(l=10, r=10, t=30, b=10)
        )
        fig = dict(data=dataplot, layout=layout)
        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_pca2djson.txt"), 'w') as outfile:
            json.dump(graphJSON, outfile)

    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_pca2djson.txt")) as json_file:
            graphJSON = json.load(json_file)

    end = time.time()
    print("create_PCA2D", end - start)
    return graphJSON


def create_PCA(filename, update):
    print("function: create_PCA()")
    start = time.time()

    graphJSON = ""

    if update == True:
        data = read_genedataset(filename)
        data = data.transpose()
        X = data.iloc[1:, :].values
        y = np.array(sampletype)

        from sklearn.preprocessing import StandardScaler
        X_std = StandardScaler().fit_transform(X)

        from sklearn.decomposition import PCA as sklearnPCA
        sklearn_pca = sklearnPCA(n_components=3)
        Y_sklearn = sklearn_pca.fit_transform(X_std)

        dataplot = []
        colors = {'Normal': '#4db6ac',
                  'Perturbation': '#ef6c00'}

        for name, col in zip(('Normal', 'Perturbation'), colors.values()):
            trace = dict(
                type='scatter3d',
                x=Y_sklearn[y == name, 0],
                y=Y_sklearn[y == name, 1],
                z=Y_sklearn[y == name, 2],
                mode='markers',
                name=name,
                marker=dict(
                    color=col,
                    size=12,
                    line=dict(
                        color='rgba(217, 217, 217, 0.14)',
                        width=0.5),
                    opacity=0.8)
            )
            dataplot.append(trace)

        layout = dict(scene=dict(
            xaxis=dict(title='PC1', showline=False),
            yaxis=dict(title='PC2', showline=False),
            zaxis=dict(title='PC3', showline=False)),
            autosize=True,
            height=500,
            margin=dict(l=10, r=10, t=30, b=10)
        )
        fig = dict(data=dataplot, layout=layout)
        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_pcajson.txt"), 'w') as outfile:
            json.dump(graphJSON, outfile)

    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_pcajson.txt")) as json_file:
            graphJSON = json.load(json_file)

    end = time.time()
    print("create_PCA", end - start)
    return graphJSON


def create_clustergrammer(filename, update):
    print("function: create_clustergrammer()")
    start = time.time()

    link = ""

    if update == True:
        df = read_genedataset(filename)
        # Construct multiple-index header
        arrays = [np.array(['Sample: ' + i for i in df.columns if i != geneheader]),
                  np.array(sampletype)]
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples)
        # Construct index
        df[geneheader] = 'Gene: ' + df[geneheader]
        # Set index and header
        df.set_index(geneheader, inplace=True)
        df.columns = index
        # Normalization
        lsize = [np.sum(df.iloc[:, i]) for i in range(len(df.columns))]
        data = (df / df.sum()) * 10 ** 6
        data = data.fillna(0)
        data = np.log2(data + (2 / (np.median(lsize) / 10 ** 6)))
        # Sorted by variance of each gene and filter only top 1500 highest variance (how distribute from mean)
        # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.var.html
        data2 = data.loc[data.var(axis=1).sort_values(ascending=False).index[:2000]]
        # Apply z-score to see how far from mean
        data3 = data2.T.apply(ss.zscore, axis=0).T
        # Create text file
        data3.to_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename + '_cluster.txt'), sep='\t', header=True,
                     index=True)
        # Upload to API and get link
        upload_url = 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'
        r = requests.post(upload_url, files={
            'file': open(os.path.join(server.config['UPLOAD_FOLDER'], filename + '_cluster.txt'), 'rb')})
        link = r.text
        print(link)

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_linkcluster.txt"), 'w') as outfile:
            outfile.write(link)

    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_linkcluster.txt")) as infile:
            link = infile.read()

    end = time.time()
    print("create_clustergrammer", end - start)
    return link


def rintegrate_test(filename, update):
    print("function: rintegrate_test()")
    start = time.time()

    if update == True:
        pandas2ri.activate()
        data = read_genedataset(filename)
        data_distinct = pandas2ri.ri2py(r.headofdataframe(filename, pandas2ri.py2ri(data), data.columns[0], sampletype))

    # print(data_distinct)
    end = time.time()
    print("rintegrate_test", end - start)


def gene_enrichment(filename, update):
    print("function: gene_enrichment()")
    start = time.time()

    graphList = []

    if update == True:
        # gene_up = [line.rstrip('\n') for line in
        #            open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_Up_genes.txt"))]
        # gene_down = [line.rstrip('\n') for line in
        #              open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_Down_genes.txt"))]
        gene_list = pd.read_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_signature.txt"))
        print(gene_list)
        gene_up = gene_list.iloc[:500, 0]
        gene_down = gene_list.iloc[-500:, 0]
        print(gene_up)
        print(gene_down)

        pathway_libraries = get_pathwaydb()
        graphList = []

        user_up = getuser(gene_up, 'Up')
        user_down = getuser(gene_down, 'Down')
        for p in pathway_libraries:
            pathway_up = enrichr(user_up, 'Up', filename, p)
            pathway_down = enrichr(user_down, 'Down', filename, p)
            graphList.append(plotenrichf(pathway_up, pathway_down))

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_enrichgraph.txt"), "w") as f:
            for s in graphList:
                f.write(str(s) + "\n")
    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_enrichgraph.txt"), "r") as f:
            for line in f:
                graphList.append(line)

    end = time.time()
    print("gene_enrichment", end - start)
    return graphList


def getuser(genelist, type):
    # Analyze gene list
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    genes_str = '\n'.join(genelist)
    # print(genes_str)
    description = 'Example gene list ' + type
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    time.sleep(1)
    user = json.loads(response.text)
    print(user)
    return user


def enrichr(useriinput, type, filename, method):
    user = useriinput
    # Get enrichment results
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'
    user_list_id = user.get('userListId')
    filename_pathway = os.path.join(server.config['UPLOAD_FOLDER'],
                                    filename + "_" + type + "_" + method + "_pathways.txt")
    gene_set_library = method

    url = ENRICHR_URL + query_string % (user_list_id, filename_pathway, gene_set_library)
    response = requests.get(url, stream=True)

    with open(filename_pathway, 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)

    # df_pathway = pd.read_csv(filename_pathway, sep="\t")
    # df_pathway = df_pathway.sort_values(by=['Combined Score'], ascending=True)
    # df_pathway = df_pathway.iloc[:10, :]
    # df_pathway = df_pathway.sort_values(by=['Combined Score'], ascending=False)
    df_pathway = pd.read_csv(filename_pathway, sep="\t")
    df_pathway = df_pathway.sort_values(by=['Combined Score'], ascending=False)
    df_pathway = df_pathway.iloc[:10, :]
    df_pathway = df_pathway.sort_values(by=['Combined Score'], ascending=True)

    return df_pathway


def plotenrichf(df_up, df_down):
    # pvlog_up = [-math.log10(x) for x in df_up['Combined Score'].tolist()]
    # pvlog_down = [-math.log10(x) for x in df_down['Combined Score'].tolist()]
    pvlog_up = [x for x in df_up['Combined Score'].tolist()]
    pvlog_down = [x for x in df_down['Combined Score'].tolist()]

    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.05, subplot_titles=("Up-regulated Pathways", "Down-regulated Pathways"))

    columns_up = df_up.columns.tolist()
    valueslist_up = df_up.values.tolist()
    hover_up = []
    for i in range(len(valueslist_up)):
        temp = []
        for j in range(len(valueslist_up[i])):
            if (columns_up[j] == 'P-value' or columns_up[j] == 'Adjusted P-value' or columns_up[j] == 'Combined Score'):
                temp.append("<b>{0}</b>: {1:.2}<br>".format(columns_up[j], valueslist_up[i][j]))
            elif (columns_up[j] == 'Genes'):
                geneoverlap = valueslist_up[i][j].split(';')
                geneoverlap = [k + "; " for k in geneoverlap]
                genebr = [x for y in (geneoverlap[i:i + 5] + ['<br>'] * (i < len(geneoverlap) - 4) for i in
                                      xrange(0, len(geneoverlap), 5)) for x in y]
                temp.append("<b>{} {}</b>: {}<br>".format(len(geneoverlap), columns_up[j], "".join(genebr)))
            elif (columns_up[j] == 'Term'):
                temp.append("<b>{}</b>: {}<br>".format(columns_up[j], valueslist_up[i][j]))
        hover_up.append(temp)
    hover_up_join = []
    for i in hover_up:
        hover_up_join.append(''.join(i))

    barfig = go.Bar(
        x=pvlog_up,
        y=df_up['Term'].tolist(),
        orientation='h',
        hovertext=hover_up_join,
        name="",
        marker=dict(color="#ff3030")
    )
    fig.append_trace(barfig, 1, 1)

    text = go.Scatter(
        x=[max(barfig['x']) / 20 for x in range(len(barfig['y']))],
        y=df_up['Term'].tolist(),
        mode='text',
        hoverinfo='none',
        showlegend=False,
        text=df_up['Term'].tolist(),
        textposition="middle right",
        textfont={'color': 'black'}
    )
    fig.append_trace(text, 1, 1)

    columns_down = df_down.columns.tolist()
    valueslist_down = df_down.values.tolist()
    hover_down = []
    for i in range(len(valueslist_down)):
        temp = []
        for j in range(len(valueslist_down[i])):
            if (columns_down[j] == 'P-value' or columns_down[j] == 'Adjusted P-value' or columns_down[
                j] == 'Combined Score'):
                temp.append("<b>{0}</b>: {1:.2}<br>".format(columns_down[j], valueslist_down[i][j]))
            elif (columns_down[j] == 'Genes'):
                geneoverlap = valueslist_down[i][j].split(';')
                geneoverlap = [k + "; " for k in geneoverlap]
                genebr = [x for y in (geneoverlap[i:i + 5] + ['<br>'] * (i < len(geneoverlap) - 4) for i in
                                      xrange(0, len(geneoverlap), 5)) for x in y]
                temp.append("<b>{} {}</b>: {}<br>".format(len(geneoverlap), columns_down[j], "".join(genebr)))
            elif (columns_down[j] == 'Term'):
                temp.append("<b>{}</b>: {}<br>".format(columns_down[j], valueslist_down[i][j]))
        hover_down.append(temp)
    hover_down_join = []
    for i in hover_down:
        hover_down_join.append(''.join(i))

    barfig2 = go.Bar(
        x=pvlog_down,
        y=df_down['Term'].tolist(),
        orientation='h',
        hovertext=hover_down_join,
        name="",
        marker=dict(color="#00bfff")
    )

    fig.append_trace(barfig2, 1, 2)

    text2 = go.Scatter(
        x=[max(barfig2['x']) / 20 for x in range(len(barfig2['y']))],
        y=df_down['Term'].tolist(),
        mode='text',
        hoverinfo='none',
        showlegend=False,
        text=df_down['Term'].tolist(),
        textposition="middle right",
        textfont={'color': 'black'}
    )
    fig.append_trace(text2, 1, 2)

    fig.update_layout(
        margin=dict(l=10, r=10, t=30, b=10),
        autosize=True,
        showlegend=False
    )
    fig.update_xaxes(title_text='Enrichment Score (Combined Score)')
    fig.update_yaxes(showticklabels=False)

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON


def pathwayinfo(filename, update):
    print("function: pathwayinfo()")
    start = time.time()

    jsonfiles = ""
    pathwaystotargetsall = pd.DataFrame()
    typelist = ['Up', 'Down']
    wordcloudlist = []

    if update == True:

        for t in typelist:
            for pl in pathway_libraries:
                pathwaystotargets = pd.DataFrame()
                df = pd.read_csv(
                    os.path.join(server.config['UPLOAD_FOLDER'], filename + "_" + t + "_" + pl + "_pathways.txt"),
                    sep="\t")

                if pl == 'KEGG_2019_Human':
                    pathwaystotargets['Pathway_Name'] = df['Term'][:10]
                    pathwaystotargets['Pathway_DB'] = pl
                    pathwaystotargets['Gene_Regulation'] = t

                    pathwayid = []
                    for p in pathwaystotargets['Pathway_Name']:
                        r = requests.get("http://rest.kegg.jp/find/pathway/" + p)
                        termid = r.text[r.text.index('path:map') + 8:r.text.index('\t')]
                        pathwayid.append("hsa" + termid)
                        time.sleep(0.5)
                    pathwaystotargets['Pathway_Id'] = pathwayid

                    uniprotralated = []
                    uniprotnumberof = []
                    for i in tqdm(pathwaystotargets['Pathway_Id']):
                        r = requests.get("http://rest.kegg.jp/link/hsa/path:" + i)
                        uniprots = []
                        geneentry_list = []
                        prefix = 'path:' + i + '\t'
                        for c in r.text.splitlines(True):
                            if (len(c) >= len(prefix)):
                                geneentry = c[c.index(prefix) + len(prefix):c.index('\n')]
                                geneentry_list.append(geneentry)
                        chunks = [geneentry_list[x:x + 50] for x in range(0, len(geneentry_list), 50)]
                        for c in chunks:
                            r2 = requests.get("http://rest.kegg.jp/conv/uniprot/" + str("+".join(c)))
                            prefix2 = '\tup:'
                            genetemp = ''
                            for c2 in r2.text.splitlines(True):
                                if (len(c2) >= len(prefix2) and c2[:c2.index('\t')] != genetemp):
                                    uniprotid = c2[c2.index(prefix2) + len(prefix2):c2.index('\n')]
                                    uniprots.append(uniprotid)
                                    genetemp = c2[:c2.index('\t')]
                                else:
                                    genetemp = ''
                            time.sleep(0.5)

                        # import multiprocessing as mp
                        # from functools import partial
                        # print("Number of processors: ", mp.cpu_count())
                        # pool = mp.Pool(mp.cpu_count())
                        # uniprotid_list = pool.map(partial(keggpathwaytotargets, i=i), r.text.splitlines(True))
                        # pool.close()
                        # pool.join()
                        # print(type(uniprotid_list))
                        # uniprots.append(uniprotid_list)

                        uniprotnumberof.append(len(list(set(uniprots))))
                        uniprotralated.append(list(set(uniprots)))
                    pathwaystotargets['Uniprot_NumberOf'] = uniprotnumberof
                    pathwaystotargets['Uniprot_Related'] = uniprotralated

                    pathwaystotargetsall = pathwaystotargetsall.append(pathwaystotargets)
                    commontargets = find_commontargets(pathwaystotargets)
                    wordcloudlist.append(commontargets)

                if pl == 'Reactome_2016':

                    pathwaystotargets['Pathway_Name'] = df['Term'][:10]
                    pathwaystotargets['Pathway_DB'] = pl
                    pathwaystotargets['Gene_Regulation'] = t

                    pathwayid = []
                    for p in pathwaystotargets['Pathway_Name']:
                        termid = p[p.index('R-HSA-'):]
                        pathwayid.append(termid)
                    pathwaystotargets['Pathway_Id'] = pathwayid

                    uniprotralated = []
                    uniprotnumberof = []
                    for i in pathwaystotargets['Pathway_Id']:
                        r = requests.get('https://reactome.org/ContentService/data/participants/' + i)
                        r_json = JSON(r.text)
                        uniprots = []
                        prefix = 'UniProt:'
                        for i in r_json.data:
                            if (len(i['refEntities']) >= 0):
                                for j in i['refEntities']:
                                    if ('UniProt' in j['displayName']):
                                        uniprotid = j['displayName'][
                                                    j['displayName'].index(prefix) + len(prefix):j['displayName'].index(
                                                        ' ')]
                                        uniprots.append(uniprotid)
                        uniprotnumberof.append(len(list(set(uniprots))))
                        uniprotralated.append(list(set(uniprots)))
                    pathwaystotargets['Uniprot_NumberOf'] = uniprotnumberof
                    pathwaystotargets['Uniprot_Related'] = uniprotralated

                    pathwaystotargetsall = pathwaystotargetsall.append(pathwaystotargets)
                    commontargets = find_commontargets(pathwaystotargets)
                    wordcloudlist.append(commontargets)

                if pl == 'WikiPathways_2019_Human':

                    pathwaystotargets['Pathway_Name'] = df['Term'][:10]
                    pathwaystotargets['Pathway_DB'] = pl
                    pathwaystotargets['Gene_Regulation'] = t

                    pathwayid = []
                    for p in pathwaystotargets['Pathway_Name']:
                        termid = p[p.index(' WP'):]
                        pathwayid.append(termid)
                    pathwaystotargets['Pathway_Id'] = pathwayid

                    import xml.etree.ElementTree as ET
                    uniprotralated = []
                    uniprotnumberof = []
                    for i in pathwaystotargets['Pathway_Id']:
                        r = requests.get('http://webservice.wikipathways.org/getXrefList?pwId=' + i + '&code=S')
                        uniprots = []
                        if r.text != '':
                            r_list = r.text.replace("ns1:", "")
                            root = ET.fromstring(r_list)
                            if root.findall('xrefs') is not None:
                                for x in root.findall('xrefs'):
                                    if x is not None:
                                        if x.text[0] == 'O' or x.text[0] == 'P' or x.text[0] == 'Q':
                                            uniprots.append(x.text)
                        uniprotnumberof.append(len(list(set(uniprots))))
                        uniprotralated.append(list(set(uniprots)))
                    pathwaystotargets['Uniprot_NumberOf'] = uniprotnumberof
                    pathwaystotargets['Uniprot_Related'] = uniprotralated

                    pathwaystotargetsall = pathwaystotargetsall.append(pathwaystotargets)
                    commontargets = find_commontargets(pathwaystotargets)
                    wordcloudlist.append(commontargets)

                print(t, pl)

        jsonfiles = json.loads(pathwaystotargetsall.to_json(orient='records'))

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_pathwayinfo.txt"), 'w') as outfile:
            json.dump(jsonfiles, outfile)

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_targetcloud.txt"), "w") as f:
            for s in wordcloudlist:
                f.write(str(s) + "\n")

    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_pathwayinfo.txt")) as json_file:
            jsonfiles = json.load(json_file)

    end = time.time()
    print("pathwayinfo", end - start)
    return jsonfiles

# def keggpathwaytotargets(c, i):
#     prefix = 'path:' + i + '\t'
#     print(c)
#     if (len(c) >= len(prefix)):
#         geneentry = c[c.index(prefix) + len(prefix):c.index('\n')]
#         r2 = requests.get("http://rest.kegg.jp/conv/uniprot/" + geneentry)
#         prefix2 = geneentry + '\tup:'
#         for c2 in r2.text.splitlines(True):
#             if (len(c2) >= len(prefix2)):
#                 uniprotid = c2[c2.index(prefix2) + len(prefix2):c2.index('\n')]
#                 break
#     time.sleep(0.5)
#     return uniprotid

def find_commontargets(pathwaystotargets):
    # Common Targets
    alluniprot = pd.DataFrame()
    uniprot = []
    for p in pathwaystotargets['Uniprot_Related']:
        for u in p:
            uniprot.append(u.strip())
    alluniprot['Uniprot'] = uniprot

    alluniprot_distinct = alluniprot.groupby(['Uniprot'], as_index=False).size().reset_index().rename(
        columns={0: 'Numbers'})
    print(alluniprot_distinct.shape)
    usort = alluniprot_distinct.sort_values(by=['Numbers'], ascending=False).reset_index(drop=True)

    allcommon = []
    for i in range(len(usort['Uniprot'])):
        common = []
        for j in range(len(pathwaystotargets['Uniprot_Related'])):
            if usort['Uniprot'][i] in pathwaystotargets['Uniprot_Related'][j]:
                common.append(pathwaystotargets['Pathway_Id'][j])
        allcommon.append(common)
    usort['Pathways'] = allcommon
    usort = usort.rename(columns={'Uniprot': 'x', 'Numbers': 'value', 'Pathways': 'pathways'})
    usort = usort.to_json(orient='records')[1:-1]
    return usort


def create_wordcloud(filename):
    cloudlist = []
    with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_targetcloud.txt"), "r") as f:
        for line in f:
            cloudlist.append(line)

    return cloudlist

def create_datatable_target(type):
    print("function: create_datatable_target()")
    start = time.time()

    naming = ['CommonPathways_', 'Pathways_', 'PathwaysDB_']
    columnnaming = [n + p for p in pathwaylibrariesloop for n in naming]

    html = ""
    html += "<table id='table_div_" + type + "' class='table display row-border'>"
    html += "<thead><tr><th class='firstth'></th>"
    html += "<th>Uniprot</th>"
    html += "<th>ProteinName</th>"
    html += "<th>Reviewed</th>"
    html += "<th>GeneName</th>"
    for c in columnnaming:
        html += "<th>" + c + "</th>"
    html += "</tr></thead><tbody></tbody><tfoot>"
    html += "<tr><th></th>"
    html += "<th>Uniprot</th>"
    html += "<th>ProteinName</th>"
    html += "<th>Reviewed</th>"
    html += "<th>GeneName</th>"
    for c in columnnaming:
        html += "<th>" + c + "</th>"
    html += "</tr></tfoot></table>"

    end = time.time()
    print("create_datatable_target", end - start)
    return html


def create_tabletarget(filename, update):
    print("function: create_tabletarget_targetvalues()")
    start = time.time()

    result_tabletarget = []

    if update == True:
        cloudlist = []
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_targetcloud.txt"), "r") as f:
            for line in f:
                cloudlist.append(line)
        # cloudlist

        print("function: create_tabletarget_targetvalues() Up")
        start = time.time()
        alluniprot_up = combine_uniprot(cloudlist, 'Up')
        end = time.time()
        print("create_tabletarget_targetvalues Up", end - start)
        alluniprot_up.to_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget_up.txt"),
                             sep='\t', encoding='utf-8')
        # alluniprot_up = pd.read_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget_up.txt"),
        #                      sep='\t', encoding='utf-8')
        # alluniprot_up = combine_uniprot_down(cloudlist, alluniprot_up, 'Up')
        print("function: create_tabletarget_targetvalues() Down")
        start = time.time()
        alluniprot_down = combine_uniprot(cloudlist, 'Down')
        end = time.time()
        print("create_tabletarget_targetvalues Down", end - start)
        alluniprot_down.to_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget_down.txt"),
                               sep='\t', encoding='utf-8')
        # alluniprot_down = pd.read_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget_down.txt"),
        #                      sep = '\t', encoding = 'utf-8')

        Row_list = []

        df = alluniprot_up
        # df = alluniprot_up.drop([alluniprot_up.columns[0]] ,  axis='columns')
        # Iterate over each row
        for index, rows in df.iterrows():
            # append the list to the final list
            Row_list.append(list(rows.values))

        result_tabletarget.append(Row_list)

        Row_list = []

        df = alluniprot_down
        # df = alluniprot_down.drop([alluniprot_down.columns[0]] ,  axis='columns')
        # Iterate over each row
        for index, rows in df.iterrows():
            # append the list to the final list
            Row_list.append(list(rows.values))

        result_tabletarget.append(Row_list)

        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget.txt"), "w") as f:
            for s in result_tabletarget:
                f.write(str(s) + "\n")

        resulttabletemp = []
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget.txt"), "r") as f:
            for line in f:
                i = 1
                temp_tabletarget = []
                newl = list(line.split('], ['))
                for l in newl:
                    ll = l.replace(" nan", " '-'")
                    ll = ll.replace("[[", "").replace("]]", "")
                    ll = "["+str(i)+"," + ll + "]"
                    temp_tabletarget.append(ast.literal_eval(ll))
                    i+=1
                resulttabletemp.append(temp_tabletarget)

        result_tabletarget = resulttabletemp

    else:
        with open(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget.txt"), "r") as f:
            for line in f:
                i = 1
                temp_tabletarget = []
                newl = list(line.split('], ['))
                for l in newl:
                    ll = l.replace(" nan", " '-'")
                    ll = ll.replace("[[", "").replace("]]", "")
                    ll = "["+str(i)+"," + ll + "]"
                    temp_tabletarget.append(ast.literal_eval(ll))
                    i+=1
                result_tabletarget.append(temp_tabletarget)

    end = time.time()
    print("create_tabletarget_targetvalues", end - start)

    return result_tabletarget


def combine_uniprot(cloudlist, type):
    # pathwaylibrariesloop = ['KEGG_2019_Human', 'WikiPathways_2019_Human', 'Reactome_2016']
    naming = ['CommonPathways_', 'Pathways_', 'PathwaysDB_']
    alluniprot = pd.DataFrame()

    for c in range(len(pathwaylibrariesloop)):
        if type == 'Up':
            data = pd.read_json(cloudlist[c], lines=True)
        elif type == 'Down':
            data = pd.read_json(cloudlist[c + len(pathwaylibrariesloop)], lines=True)
        data['pathwaysdb'] = pathwaylibrariesloop[c]
        if c == 0:
            alluniprot = data
            alluniprot.columns = ['x'] + [str(col) + '_' + pathwaylibrariesloop[c] for col in alluniprot.columns if
                                          col != 'x']
        else:
            alluniprot = pd.merge(alluniprot, data, on='x', how='outer',
                                  suffixes=('_' + pathwaylibrariesloop[c - 1], '_' + pathwaylibrariesloop[c]))

    columnnaming = [n + p for p in pathwaylibrariesloop for n in naming]
    alluniprot.columns = ['Uniprot'] + columnnaming

    # urluniprot = '{http://uniprot.org/uniprot}'
    # reviewed = []
    # unigene = []
    # uniname = []
    # i = 1
    # for u in alluniprot['Uniprot']:
    #     uniprotxml = requests.get("https://www.uniprot.org/uniprot/" + u.strip() + ".xml")
    #     if uniprotxml.text != '':
    #         root = ET.fromstring(uniprotxml.text)
    #         rv = root[0].attrib['dataset']
    #         gn = ''
    #         if root.find(urluniprot + 'entry').find(urluniprot + 'gene') is not None:
    #             gn = root.find(urluniprot + 'entry').find(urluniprot + 'gene').find(urluniprot + 'name').text
    #         else:
    #             gn = 'No Gene'
    #         fn = ''
    #         if rv == 'Swiss-Prot':
    #             fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
    #                 urluniprot + 'recommendedName').find(
    #                 urluniprot + 'fullName').text
    #         elif rv == 'TrEMBL':
    #             if root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
    #                     urluniprot + 'submittedName') is not None:
    #                 fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
    #                     urluniprot + 'submittedName').find(urluniprot + 'fullName').text
    #             else:
    #                 fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein')[0].find(
    #                     urluniprot + 'fullName').text
    #     else:
    #         rv = "Obsolete Entry"
    #         gn = "Obsolete Entry"
    #         fn = "Obsolete Entry"
    #     reviewed.append(rv)
    #     unigene.append(gn)
    #     uniname.append(fn)
    #     print(u, i)
    #     i += 1
    # alluniprot['Reviewed'] = reviewed
    # alluniprot['GeneName'] = unigene
    # alluniprot['ProteinName'] = uniname
    # alluniprot

    finaldataframe = pd.DataFrame(columns=['Uniprot', 'Reviewed', 'GeneName', 'ProteinName'])
    # keys = []
    # count = 0
    chunks = [alluniprot['Uniprot'][x:x + 100] for x in xrange(0, len(alluniprot['Uniprot']), 100)]

    import multiprocessing as mp
    print("Number of processors: ", mp.cpu_count())
    pool = mp.Pool(mp.cpu_count())
    dfframe = pool.map(requestuniprot, chunks)

    # for i in tqdm(alluniprot['Uniprot']):
    #     keys.append(i.strip())
    #     count += 1
    #     if (count % 100) == 0 or count == len(alluniprot['Uniprot']):
    #         requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?accession=" + ",".join(keys)
    #         r = requests.get(requestURL, headers={"Accept": "application/json"})
    #         keys = []
    #         df = json_normalize(json.loads(r.text))
    #         # print(df.columns)
    #         u = df['accession']
    #         rv = df['info.type']
    #         fn = pd.Series([i[0]['name']['value'] if str(i) != "nan" else i for i in df['gene']])
    #         gn = df['protein.recommendedName.fullName.value']
    #         frame = {'Uniprot': u, 'Reviewed': rv, 'ProteinName': fn, 'GeneName': gn}
    #         dfframe = pd.DataFrame(frame)
    #         finaldataframe = finaldataframe.append(dfframe)

    pool.close()
    for d in dfframe:
        finaldataframe = finaldataframe.append(d)

    alluniprot = pd.merge(alluniprot, finaldataframe, how='outer', left_on='Uniprot', right_on='Uniprot')

    # Using Multiprocessing Techniques
    # import multiprocessing as mp
    # print("Number of processors: ", mp.cpu_count())
    # pool = mp.Pool(mp.cpu_count())
    # results = pool.map(requestuniprot, alluniprot['Uniprot'])
    # pool.close()
    # dfresults = pd.DataFrame(results, columns=['Uniprot', 'Reviewed', 'GeneName', 'ProteinName'])
    # alluniprot = pd.merge(alluniprot, dfresults, how='inner', left_on='Uniprot', right_on='Uniprot')

    alluniprot = alluniprot[['Uniprot', 'ProteinName', 'Reviewed', 'GeneName'] + columnnaming]

    return alluniprot

# def requestuniprot(u):
#     urluniprot = '{http://uniprot.org/uniprot}'
#     uniprotxml = requests.get("https://www.uniprot.org/uniprot/" + u.strip() + ".xml")
#     if uniprotxml.text != '':
#         root = ET.fromstring(uniprotxml.text)
#         rv = root[0].attrib['dataset']
#         gn = ''
#         if root.find(urluniprot + 'entry').find(urluniprot + 'gene') is not None:
#             gn = root.find(urluniprot + 'entry').find(urluniprot + 'gene').find(urluniprot + 'name').text
#         else:
#             gn = 'No Gene'
#         fn = ''
#         if rv == 'Swiss-Prot':
#             fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
#                 urluniprot + 'recommendedName').find(
#                 urluniprot + 'fullName').text
#         elif rv == 'TrEMBL':
#             if root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
#                     urluniprot + 'submittedName') is not None:
#                 fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
#                     urluniprot + 'submittedName').find(urluniprot + 'fullName').text
#             else:
#                 fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein')[0].find(
#                     urluniprot + 'fullName').text
#     else:
#         rv = "Obsolete Entry"
#         gn = "Obsolete Entry"
#         fn = "Obsolete Entry"
#     # reviewed.append(rv)
#     # unigene.append(gn)
#     # uniname.append(fn)
#     print(u)
#     return u, rv, gn, fn

def requestuniprot(chunks):
    # print(chunks)
    keys = [x.strip(' ') for x in chunks]
    keys = [x for x in keys if len(x) > 5]
    # print(",".join(keys))
    requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?accession=" + ",".join(keys)
    r = requests.get(requestURL, headers={"Accept": "application/json"})
    keys = []
    df = json_normalize(json.loads(r.text))
    u = df['accession']
    rv = df['info.type']
    gn = pd.Series([i[0]['name']['value'] if str(i) != "nan" and 'name' in i[0].keys() else np.nan for i in df['gene']])
    fn = df['protein.recommendedName.fullName.value']
    frame = {'Uniprot': u, 'Reviewed': rv, 'ProteinName': fn, 'GeneName': gn}
    dfframe = pd.DataFrame(frame)
    return dfframe


# def combine_uniprot_down(cloudlist, alluniprot_up, type):
#     # pathwaylibrariesloop = ['KEGG_2019_Human', 'WikiPathways_2019_Human', 'Reactome_2016']
#     naming = ['CommonPathways_', 'Pathways_', 'PathwaysDB_']
#     alluniprot = pd.DataFrame()
#
#     for c in range(len(pathwaylibrariesloop)):
#         if type == 'Up':
#             data = pd.read_json(cloudlist[c], lines=True)
#         elif type == 'Down':
#             data = pd.read_json(cloudlist[c + len(pathwaylibrariesloop)], lines=True)
#         data['pathwaysdb'] = pathwaylibrariesloop[c]
#         if c == 0:
#             alluniprot = data
#             alluniprot.columns = ['x'] + [str(col) + '_' + pathwaylibrariesloop[c] for col in alluniprot.columns if
#                                           col != 'x']
#         else:
#             alluniprot = pd.merge(alluniprot, data, on='x', how='outer',
#                                   suffixes=('_' + pathwaylibrariesloop[c - 1], '_' + pathwaylibrariesloop[c]))
#
#     columnnaming = [n + p for p in pathwaylibrariesloop for n in naming]
#     alluniprot.columns = ['Uniprot'] + columnnaming

    # urluniprot = '{http://uniprot.org/uniprot}'
    # reviewed = []
    # unigene = []
    # uniname = []
    # i = 1
    # for u in alluniprot['Uniprot']:
    #     if u in list(alluniprot_up['Uniprot']):
    #         print(i, u)
    #         row = alluniprot_up.loc[alluniprot_up['Uniprot'] == u].index[0]
    #         rv = alluniprot_up.loc[row, 'Reviewed']
    #         gn = alluniprot_up.loc[row, 'GeneName']
    #         fn = alluniprot_up.loc[row, 'ProteinName']
    #     else:
    #         uniprotxml = requests.get("https://www.uniprot.org/uniprot/" + u.strip() + ".xml")
    #         print(u)
    #         if uniprotxml.text != '' and uniprotxml is not None:
    #             root = ET.fromstring(uniprotxml.text)
    #             rv = root[0].attrib['dataset']
    #             gn = ''
    #             if root.find(urluniprot + 'entry').find(urluniprot + 'gene') is not None:
    #                 gn = root.find(urluniprot + 'entry').find(urluniprot + 'gene').find(urluniprot + 'name').text
    #             else:
    #                 gn = 'No Gene'
    #             fn = ''
    #             if rv == 'Swiss-Prot':
    #                 fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
    #                     urluniprot + 'recommendedName').find(
    #                     urluniprot + 'fullName').text
    #             elif rv == 'TrEMBL':
    #                 if root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
    #                         urluniprot + 'submittedName') is not None:
    #                     fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein').find(
    #                         urluniprot + 'submittedName').find(urluniprot + 'fullName').text
    #                 else:
    #                     fn = root.find(urluniprot + 'entry').find(urluniprot + 'protein')[0].find(
    #                         urluniprot + 'fullName').text
    #         else:
    #             rv = "Obsolete Entry"
    #             gn = "Obsolete Entry"
    #             fn = "Obsolete Entry"
    #     reviewed.append(rv)
    #     unigene.append(gn)
    #     uniname.append(fn)
    #     print(u, i)
    #     i += 1
    # alluniprot['Reviewed'] = reviewed
    # alluniprot['GeneName'] = unigene
    # alluniprot['ProteinName'] = uniname
    # alluniprot



    # alluniprot = alluniprot[['Uniprot', 'ProteinName', 'Reviewed', 'GeneName'] + columnnaming]
    #
    # return alluniprot


def dict_target(filename):
    alluniprot_up = pd.read_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget_up.txt"),
                         sep='\t', encoding='utf-8')
    alluniprot_down = pd.read_csv(os.path.join(server.config['UPLOAD_FOLDER'], filename + "_datavaluestarget_down.txt"),
                           sep='\t', encoding='utf-8')

    alluniprot = alluniprot_up.append(alluniprot_down, ignore_index=True)
    alluniprot = alluniprot.drop(columns=['Unnamed: 0'])
    alluniprot_json = json.loads(alluniprot.to_json(orient='records'))
    return alluniprot_json