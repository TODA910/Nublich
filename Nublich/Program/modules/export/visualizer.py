import plotly
import plotly.express as px
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter, gaussian_filter1d

# prepare dataframe for summary table
def get_sum_df(df):
    num_barcode = df['barcode'].value_counts().shape[0]
    reads = len(df)
    bases = df['read_length'].sum()
    mean_phred_score = round(sum([float(i) for i in df['mean_phred_score']])/reads, 2)
    mean_read_length = round(bases/reads, 2)

    dic = {'Reads':reads,
           'Bases':bases,
           'Mean read length':mean_read_length,
           'Mean phred score':mean_phred_score,
           'Barcodes':num_barcode}
    sum_df = pd.DataFrame(dic, index=['row1'])
    return sum_df

# prepare dataframe for barcode summary table
def get_barcodesum_df(df):
    all_barcode = pd.unique(df['barcode'])
    reads = []
    bases = []
    mean_read_length = []
    mean_phred_score = []
    
    grouped = df.groupby('barcode')
    
    for barcode in all_barcode:
        df_by_barcode = grouped.get_group(barcode)
        n = len(df_by_barcode)
        b = sum([float(i) for i in df_by_barcode['read_length']])
        reads.append(n)
        bases.append(b)
        mean_read_length.append(round(b/n, 2))
        mean_phred_score.append(round(sum([float(i) for i in df_by_barcode['mean_phred_score']])/n, 2))
        
    dic = {'Barcode':all_barcode,
           'Reads':reads,
           'Bases':bases,
           'Mean read length':mean_read_length,
           'Mean phred score':mean_phred_score
           }
    
    barcodesum_df = pd.DataFrame(dic)
    return barcodesum_df

# create table
def getTable(df, Table_name):
    d = go.Table(
#     columnwidth=[2,3,4,6], 
            header=dict(
                values=list(df.keys()),
                font=dict(size=15),
#                 align="left"
            ),
            cells=dict(
                values=[df[k].tolist() for k in df.keys()],
                font=dict(size=15),
                height=40
#             align = "left"
            )
        )
    fig = go.Figure(data=d)
    fig.update_layout(height = 350)
    divfig = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
#     fig.show(df)
    return {Table_name:divfig}

# create length density graph
def graph_length_density(df, graph_name):
    allrecord_length = [float(i) for i in list(df['read_length'])]
    stat = np.percentile (allrecord_length, [10,25,50,75,90])
    ## เปลี่ยนตรง data ด้วย nbins ไว้เปลี่ยนตัวแปลตอนทำฟังชั่น   อันนี้ดูด้วยว่ามันต่างกันยังไงว่า log กับ ไม่ log                       
    mindata = np.nanmin(allrecord_length)
    maxdata = np.nanmax(allrecord_length)
    # hisx = np.histogram (a=allrecord_length, bins= np.linspace (mindata,maxdata , 100))
    hisx = np.histogram (a=allrecord_length, bins=np.logspace (np.log10(mindata), np.log10(maxdata)+0.1, 200))
    
    x = [hisx[1],[stat[0],stat[0]], [stat[1],stat[1]], [stat[2],stat[2]], [stat[3],stat[3]], [stat[4],stat[4]]]

    y_count = gaussian_filter1d(hisx[0],sigma=2)
    y_max = y_count.max()
    y = [y_count,[0,y_max], [0,y_max], [0,y_max], [0,y_max], [0,y_max]]

    name = ["Density", "10%", "25%", "Median", "75%", "90%"]
    text = ["",
                ["", "10%<br>{:,.2f}".format(stat[0])],
                ["", "25%<br>{:,.2f}".format(stat[1])],
                ["", "Median<br>{:,.2f}".format(stat[2])],
                ["", "75%<br>{:,.2f}".format(stat[3])],
                ["", "90%<br>{:,.2f}".format(stat[4])]]
    
    common = {"mode": "lines+text",
              "hoverinfo": "skip",
              "textposition": 'top center',
              "line":  {'color':'gray','width':1,'dash': 'dot'}} 
    
    data = [go.Scatter (x=x[0], y=y[0],name=name[0], fill='tozeroy', fillcolor="red", mode='none', showlegend=True),
            go.Scatter (x=x[1], y=y[1], name=name[1], text=text[1], **common),
            go.Scatter (x=x[2], y=y[2], name=name[2], text=text[2], **common),
            go.Scatter (x=x[3], y=y[3], name=name[3], text=text[3], **common),
            go.Scatter (x=x[4], y=y[4], name=name[4], text=text[4], **common),
            go.Scatter (x=x[5], y=y[5], name=name[5], text=text[5], **common)]
    
    layout = go.Layout (hovermode = "closest",
                        plot_bgcolor="whitesmoke",
                        legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
                        width = None,
                        height = 500,
                        title = {"xref":"paper" ,"x":0.5, "xanchor":"center"},
                        xaxis = {"title":"length","type":"log", "zeroline":False, "showline":True},
                        yaxis = {"title":"Read density", "zeroline":False, "showline":True, "fixedrange":True,"range":[0, y_max+y_max/6]})
    
    fig = go.Figure (data=data, layout=layout)
    
    divfig = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
    return {graph_name:divfig}

# create phred density graph
def graph_phred_density(df, graph_name, stat_info):
    meanphred = [float(i) for i in list(df['mean_phred_score'])]
    stat = stat_info[1:6]
    mindata = stat_info[0]
    maxdata = stat_info[6]
    hisx = np.histogram (a=meanphred, bins= np.linspace (mindata,maxdata , 100))

    x = [hisx[1],[stat[0],stat[0]], [stat[1],stat[1]], [stat[2],stat[2]], [stat[3],stat[3]], [stat[4],stat[4]]]

    y_count = gaussian_filter1d(hisx[0],sigma=2)
    y_max = y_count.max()
    y = [y_count,[0,y_max], [0,y_max], [0,y_max], [0,y_max], [0,y_max]]

    name = ["Density", "10%", "25%", "Median", "75%", "90%"]
    text = ["",
                ["", "10%<br>{:,.2f}".format(stat[0])],
                ["", "25%<br>{:,.2f}".format(stat[1])],
                ["", "Median<br>{:,.2f}".format(stat[2])],
                ["", "75%<br>{:,.2f}".format(stat[3])],
                ["", "90%<br>{:,.2f}".format(stat[4])]]

    common = {"mode": "lines+text",
              "hoverinfo": "skip",
              "textposition": 'top center',
              "line":  {'color':'gray','width':1,'dash': 'dot'}}                                                                           

    data = [go.Scatter (x=x[0], y=y[0],name=name[0], fill='tozeroy', fillcolor="dodgerblue", mode='none', showlegend=True),
            go.Scatter (x=x[1], y=y[1], name=name[1], text=text[1], **common),
            go.Scatter (x=x[2], y=y[2], name=name[2], text=text[2], **common),
            go.Scatter (x=x[3], y=y[3], name=name[3], text=text[3], **common),
            go.Scatter (x=x[4], y=y[4], name=name[4], text=text[4], **common),
            go.Scatter (x=x[5], y=y[5], name=name[5], text=text[5], **common)]

    layout = go.Layout (
        hovermode = "closest",
        plot_bgcolor="whitesmoke",
        legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
        width = None,
        height = 500,
        title = {"xref":"paper" ,"x":0.5, "xanchor":"center"},
        xaxis = {"title":"Phred score","type":"linear", "zeroline":False, "showline":True},
        yaxis = {"title":"Density", "zeroline":False, "showline":True, "fixedrange":True,"range":[0, y_max+y_max/6]})

    fig = go.Figure (data=data, layout=layout)
    
    divfig = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
    return {graph_name:divfig}

# create pie chart
def graph_pie(df, graph_name):
    df = df['barcode']
    dg = df.value_counts()
    dh = dg.to_dict()

    all_barcodes = []
    num_barcode = []
    for key,value in dh.items():
        all_barcodes.append(key)
        num_barcode.append(value)
    data = go.Pie (
            labels=all_barcodes,
            values=num_barcode,
            sort=False,
            marker={"colors":["#f44f39","#fc8161","#fcaf94","#828282"]},
            name="Pie plot",
            textinfo='label+percent')
    
    fig = go.Figure (data=data)
    
    divfig = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
    return {graph_name:divfig}

def graph_length_phred(df, graph_name):
    meanphred = [float(i) for i in df['mean_phred_score'].values]
    allrecord_length = [float(i) for i in df['read_length'].values]
    
    mindata_length = np.nanmin(allrecord_length)
    maxdata_length = np.nanmax(allrecord_length)
    mindata_phred = np.nanmin(meanphred)
    maxdata_phred = np.nanmax(meanphred)
    
    binx=np.logspace (np.log10(mindata_length), np.log10(maxdata_length)+0.1, 200)
    biny=np.linspace (mindata_phred,maxdata_phred , 200)
    
    hist = np.histogram2d (x=meanphred, y=allrecord_length, bins=[biny,binx])

    z_min = np.percentile (hist[0], ( 0))
    z_max = np.percentile (hist[0], ( 100))     
    z = gaussian_filter(hist[0], sigma=2)
    data = [go.Contour (x=hist[2], y=hist[1], z=z,
                        name="Density",
                        hoverinfo="name+x+y",
                        colorscale="hot",
                        showlegend=True,
                        connectgaps=True,
                        line={"width":0})]

    layout = go.Layout (hovermode = "closest",
                        plot_bgcolor="whitesmoke",
                        legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
#                         width = 500,
                        height = 500,
                        title = {"xref":"paper" ,"x":0.5, "xanchor":"center"},
                        xaxis = {"title":"length", "showgrid":True, "zeroline":False, "showline":True, "type":"log"},
                        yaxis = {"title":"Phred score", "showgrid":True, "zeroline":False,"showline":True, "type":"linear"})

    fig = go.Figure (data=data, layout=layout)
    
    divfig = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
    return {graph_name:divfig}