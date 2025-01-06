import pandas as pd
import requests, os
from flask import request, Response
from time import sleep
from rdkit import Chem
from rdkit.Chem import PandasTools

from rdkit.Chem.Draw import rdMolDraw2D

from rdkit.Chem import Draw
from json import dumps
from flask import Flask, render_template
import math
import base64
import uuid
import hashlib
import helpers

app = Flask(__name__)

#filename = ("data/wikipedia.csv")
#filename = ("data/chembl_sars.sdf")

@app.route('/')
def index():
    list_string = ""
    file_list = os.listdir('data/')
    for file in file_list:
        list_string += f"<option value='{file}'>{file}</option>"
    return render_template('index.html', list_string=list_string)


@app.route('/table/<filename>')
def table(filename):
    filename = "data/" + filename
    uid = uuid.uuid4().hex
    if filename.endswith('.csv'):
        print('Loading csv')
        df = pd.read_csv(filename, nrows= 999999)
        #PandasTools.AddMoleculeColumnToFrame(df, 'SMILES', 'Molecule')
        df['Molecule'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x, sanitize=True) )
        # if Molecule is None then Molecule become empty molecule
        df['Molecule'] = df['Molecule'].apply(lambda Molecule: Molecule if Molecule else Chem.MolFromSmiles('') )
    elif filename.endswith('.sdf'):
        print('Loading sdf')
        #df = PandasTools.LoadSDF(filename, molColName='Molecule')
        mol_supplier = Chem.MultithreadedSDMolSupplier(filename, numWriterThreads=5)
        molecules = []
        for mol in mol_supplier:
            if mol is not None:
                props = mol.GetPropsAsDict()
                props["Title"] = mol.GetProp("_Name")
                props["Molecule"] = mol
                molecules.append(props)
        df = pd.DataFrame(molecules)

    #if("Name" not in df.columns):
    #    df['Name'] = df['Molecule'].apply(lambda x: Chem.InchiToInchiKey(Chem.MolToInchi(x, options='-KET -15T')))

    df.insert(0, 'Molecule', df.pop('Molecule'))
    
    df.reset_index(inplace=True)
    globals()[uid] = df 

    table = "<table class='maintable'>"
    for index, row in df.iterrows():
        if (index == 0 ):
            table += "<tr>"
            for column in df.columns:
                if column == 'index': 
                    table += "<th  class='stickyHeaderIndex'>&nbsp;</th>"        
                elif column == 'Molecule':
                    table += "<th  class='stickyHeaderMolecule'>Molecule</th>"     
                elif column == 'Name':
                    table += "<th  class='stickyHeaderName'>Name</th>"  
                else:
                    class_sticky = ""
                    column_name = str(column).replace('_',' ').replace('*',' ').replace('.','<br>')
                    table += "<th name='"+ str(column) + "' class='sortableHeader' onclick='sortHeader(this)'  title='Click to sort by this column'>" + column_name + "</th>"        
            table += "</tr>"    
        table += f"<tr id='{str(index)}' class='observable' ></tr>"
    table += "</table>"
    return render_template('main.html', table=table, uid = uid)
    


@app.route("/readrow/" , methods = ['POST'])
def readrow():
    session = str(request.form.get('session'))  
    index = int(str(request.form.get('index')))
    try:
        df = globals()[session]
        row = df.iloc[index] #.query('index == ' + str(index) )
        return renderRow(row, df.columns)
    except:
        return Response("", status=204)
    


def renderRow(row, columns):
    table = "<tr>"
    for column in columns:
        if ( column == "Molecule"):
            try:
                drawer = rdMolDraw2D.MolDraw2DCairo(200, 100)
                drawer.drawOptions().useBWAtomPalette()
                drawer.SetLineWidth(0.9)
                drawer.drawOptions().minFontSize = 9
                color=(245/255.0,245/255.0, 245/255.0, 245/255.0)
                drawer.drawOptions().setBackgroundColour(color)



                rdMolDraw2D.PrepareAndDrawMolecule(drawer, row['Molecule'])                
                table += f"""<td class='stickyMolecule'><div style=' height: 100px; width: 200px;'><img src='data:image/png;base64, {base64.b64encode(drawer.GetDrawingText()).decode('utf8')}' /></div></td>"""
            except Exception as ex :
                print("An exception occurred ")
                table += f"<td class='stickyMolecule'><div style=' display: table-cell;  border: 0px solid black; height: 104px; width: 204px;'><img src='data:image/png;base64,xxx' /></div></td>"
        elif ( column == "SMILES"):
            table += "<td class='smiles' >" + str(row[column]) + "</td>"
        elif ( column == "index"):
            table += "<td class='stickyIndex' >" + str(row[column]) + "</td>"
        elif ( column == "Name"):
            table += "<td class='stickyName'><a target='_blank' href='https://en.wikipedia.org/wiki/" + str(row[column]) + "'>" + str(row[column]) + "</a></td>"
        
        else:
            
            if isinstance(row[column],float):
                if math.isnan(row[column]):
                    tdvalue = ""
                #elif row[column] == 0:
                #    tdvalue = "0"
                else:
                    try:
                        tdvalue = "<span title='"+str(row[column])+"'>" + str(helpers.round_to_one_significant_decimal(row[column])) + "</span>"
                    except Exception as ex :
                        print("An exception occurred ")
                        tdvalue = "<span style='color:red' >" + str(row[column]) + "</span>"

            else:
                tdvalue = str(row[column] )
            table += "<td>" + tdvalue + "</td>"
    table += "</tr>"
    return table


