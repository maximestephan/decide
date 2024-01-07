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


app = Flask(__name__)


@app.route('/')
def search():
    return render_template('main.html')



