import pandas as pd
import requests, os
from flask import request
from time import sleep
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdFMCS 
from rdkit import DataStructs
from rdkit import RDConfig
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem import rdmolops
from rdkit.Chem import Draw
from json import dumps
from flask import Flask, render_template
import base64
import hashlib
import helpers
import matplotlib.pyplot as plt

app = Flask(__name__)


@app.route('/')
def search():
    data = [(50,100),(12.5,90),(3.125,100),(0.7813,70),(0.1953,50),(0.0488,20)]

    base64img = helpers.display_biological_area( data )
    html = f"<img  src='data:image/jpeg;base64,{base64img}' >"
    return html
    return render_template('main.html')
    "*/width='100' height='100'/*"


