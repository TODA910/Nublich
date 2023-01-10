import plotly
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from pathlib import Path
from dominate.tags import *
from dominate import document

# get "dic" from function getContent()
def createSideNav(dic):
    sideNav = div(id="mySidenav", cls='sidenav', onclick='closeNav()')
    with sideNav:
        a('Home', href='#')
        for i in dic:
            a(dic[i].keys(), href='#%s'% i)
    return sideNav.render()

# get "dic" from function getContent()
def createMain(dic):
    main = div(cls='main')
    with main:
        for i in dic:
            branch = div(cls='sec',id=i)
            with branch:
                h2(dic[i].keys(),cls='title')
                div(dic[i].values())
    return main.render().replace('&quot;','"').replace('&lt;','<').replace('&gt;','>')

def getTemplate():
    template = Path('template.html').read_text()
    return template

# get "content" from function getContent()
def fullHtml(content):
    template = getTemplate()
    sideNav = createSideNav(content)
    main = createMain(content)
    fullhtml = template.replace('forSideNav', sideNav).replace('forMain', main)
    return fullhtml

# get "fullhtml" from function fullHtml()
def writeHtml(fullhtml,filename):
    with open(f'{filename}_report.html', 'w') as report:
        report.write(fullhtml)

# "argument" must be in this form --> {'description':plotly_object}
# "content" will be in this form --> {id:{'description':plotly_object}}
def getContent(tupl):
    content = {}
    for i,arg in enumerate(tupl):
        ID = 'content'+str(i+1)
        content[ID] = arg
    return content

# "argument" must be in this form --> {'description':plotly_object}
def createReport(*argument,filename):
    content = getContent(argument)
    fullhtml = fullHtml(content)
    writeHtml(fullhtml,filename)