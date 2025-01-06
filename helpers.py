import io
import base64
import cmath
from math import log10, floor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools


def round_to_one_significant_decimal(x):
    if "e" in str(f"{x:.3g}"):
        if x > 1 or x < -1:
            return int(x)
        else:
            decimals = -floor(log10(abs(x)))
            as_str = f'{{:.{decimals}f}}'.format(x)
            chunks = []
            start_chunk = as_str.find('.') + 4
            chunks.append(as_str[:start_chunk])
            for i in range(start_chunk, len(as_str), 3):
                chunks.append(as_str[i:i+3])
            return ' '.join(chunks)

    return f"{x:.3g}"
    #return round(x, -int(floor(log10(abs(x)))))


def calculate_pBP8x(list_tupple_rate_efficacy, x):
    pBP8x = 0
    for tupple_rate_efficay in list_tupple_rate_efficacy:
        rate = tupple_rate_efficay[0]
        efficacy = tupple_rate_efficay[1]
        if efficacy >= x:
            pBP8x = -cmath.log10(rate/10 /1000 ).real 
    return pBP8x



def calculate_biological_score(list_tupple_rate_efficacy):
    pBPx = 0
    sum_threslhold = 0
    for x in range(10):
        threslhold = (x+1)*10
        pBPx += threslhold * calculate_pBP8x(list_tupple_rate_efficacy,threslhold )
        sum_threslhold += threslhold
    score = round( pBPx / sum_threslhold , 2)
    return score

def Area(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area


def display_biological_area(list_tupple_rate_efficacy):
    n = len(list_tupple_rate_efficacy) 
    xlist = []
    ylist = []
    
    for i in range(n):
        xlist.append( list_tupple_rate_efficacy[i][0] )
        ylist.append( list_tupple_rate_efficacy[i][1])

    xpolygon = xlist
    xpolygon.insert(0, list_tupple_rate_efficacy[0][0])
    xpolygon.append(list_tupple_rate_efficacy[n-1][0] )
    
    ypolygon = ylist
    ypolygon.insert(0, 0)
    ypolygon.append(0)
    
    plt.xscale('log')
    plt.xticks(ticks= xlist, labels=xlist, minor=False, rotation='vertical')
    
    #plt.xlabel('Concentration')
    plt.xlabel('')
    
    plt.gca().get_xaxis().set_inverted("inverse")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlim((max(xlist),min(xlist)))
    
    plt.yticks(ticks= ylist, label = ylist)
    #plt.ylabel('Efficacy')
    plt.ylim((0,105))
    plt.plot(xlist, ylist, marker = 'o')
    
    plt.fill(xpolygon, ypolygon, alpha=0.2, facecolor='yellow')
    #plt.figure(figsize=(2, 2))
    
    my_stringIObytes = io.BytesIO()
    plt.savefig(my_stringIObytes, format='png' )
    my_stringIObytes.seek(0)
    img_as_base64 = base64.b64encode(my_stringIObytes.read()).decode()
    plt.close()

    return img_as_base64
    

def calculate_biological_area(list_tupple_rate_efficacy):
    n = len(list_tupple_rate_efficacy)-1 
    area = 0.0
    for i in range(n):
        width = list_tupple_rate_efficacy[i][0] - list_tupple_rate_efficacy[i+1][0]
        heigth = min(list_tupple_rate_efficacy[i][1] , list_tupple_rate_efficacy[i+1][1])
        area_rectangle = width * heigth
        heigth = abs(list_tupple_rate_efficacy[i][1] - list_tupple_rate_efficacy[i+1][1])
        area_triangle  = width * heigth / 2
        area += area_rectangle + area_triangle
    return round(area / 100 ,2 ) 
    
def using_display_biological_area():
    results = None #getResults()
    html = ""
    for result in results:
    
        base64img = display_biological_area( result[1] )
        html += f"<img  src='data:image/jpeg;base64,{base64img}' >"
    return html


def cleanFile(filename):
    df = pd.read_csv(filename)

    #PandasTools.AddMoleculeColumnToFrame(df, 'SMILES', 'Molecule')
    df['Molecule'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x, sanitize=True) )
    
    # if Molecule is None then Molecule become empty molecule
    df['Molecule'] = df['Molecule'].apply(lambda Molecule: Molecule if Molecule else Chem.MolFromSmiles('') )
    #df['InchiKey'] = df['Molecule'].apply(lambda x: Chem.InchiToInchiKey(Chem.MolToInchi(x, options='-KET -15T')))
    #df = df[['InchiKey','Name','SMILES']]


    df['logP'] = df['Molecule'].apply(lambda mol: round(Descriptors.MolLogP(mol),1) )

    descrs = [Descriptors.CalcMolDescriptors(mol) for mol in df['Molecule'].values.tolist()]
    df = df.join(pd.DataFrame(descrs))
     
    df.drop(columns=['Molecule'], inplace=True)
    df.drop(columns=df.columns[0], axis=1, inplace=True)
    df.to_csv(filename + '_out.csv',index=False)
            


#cleanFile('data/wikipedia.csv')

