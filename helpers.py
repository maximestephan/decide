import io
import base64
import cmath
#from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem

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


def calculate_biological_area2(list_tupple_rate_efficacy):
    copy_list = list_tupple_rate_efficacy
    origin_list = []
    for index, tupple_rate_efficay in enumerate(copy_list):
        copy_list[index] = (  tupple_rate_efficay[0] /10 , tupple_rate_efficay[1])
        origin_list.append( (tupple_rate_efficay[0] /10 , 0))
    origin_list.reverse()
    t = tuple(copy_list + origin_list)
    polygon = Polygon(t)
    plt.plot(*polygon.exterior.xy)
    plt.show()
    return round( polygon.area / 100  ,2 )

def transform(x):
    #return x
    return cmath.log(x).real

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
    

def cleanFile(filename):
    
    df = pd.read_csv(filename)

    #PandasTools.AddMoleculeColumnToFrame(df, 'SMILES', 'Molecule')
    df['Molecule'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x, sanitize=False) )

    # if Molecule is None then Molecule become empty molecule
    df['Molecule'] = df['Molecule'].apply(lambda Molecule: Molecule if Molecule else Chem.MolFromSmiles('') )
    df['InchiKey'] = df['Molecule'].apply(lambda x: Chem.InchiToInchiKey(Chem.MolToInchi(x, options='-KET -15T')))
    df = df[['InchiKey','Name','SMILES']]

    #df.drop(columns=['Molecule'], inplace=True)
    #filtered_df = df[df['InchiKey'].isnull()]

    df.to_csv(filename + '_out.csv')
            

    #print(filtered_df)

#list_tupple_rate_efficacy = [(100, 90),(50, 20),(25, 20),(13, 0),(6, 0),(3, 0)]
#score = calculate_biological_score(list_tupple_rate_efficacy)
#area = calculate_biological_area(list_tupple_rate_efficacy)
#print(score,area)

#list_tupple_rate_efficacy = [(100, 100),(50, 100),(25, 100),(13, 100),(6, 100),(3, 100)]
#score = calculate_biological_score(list_tupple_rate_efficacy)
#area = calculate_biological_area(list_tupple_rate_efficacy)
#print(score,area)


results = [ 
(1.26,[(100,96),(50,96),(25,83)] ) , 
(1.05,[(100,73)]),
(0.78,[(100,77),(50,70),(25,62)]),
(1.50,[(100,100),(50,100),(25,88)]),
(1.44,[(100,100),(50,97),(25,84)]),
(1.21,[(100,97),(50,88),(25,84)]),
(1.43,[(100,96),(50,54),(25,61),(13,72),(6,79)]),
(3.50,[(100,100),(50,100),(25,100),(13,100),(6,100),(3,100)]),
(3.50,[(100,100),(50,100),(25,100),(13,100),(6,100),(3,100)]),
(3.50,[(100,100),(50,100),(25,100),(13,100),(6,100),(3,100)]),
(1.67,[(100,90),(50,20),(25,20),(13,0),(6,0),(3,0)]),
(0.00,[(100,0),(50,0),(25,0),(13,0),(6,0),(3,0)]),
(3.50,[(100,100),(50,100),(25,100),(13,100),(6,100),(3,100)]),
(3.50,[(100,100),(50,100),(25,100),(13,100),(6,100),(3,100)]),
(2.47,[(100,100),(50,100),(25,70),(13,20),(6,0),(3,0)]),
(3.45,[(100,100),(50,100),(25,100),(13,100),(6,100),(3,90)]),
(2.75,[(100,100),(50,100),(25,100),(13,70),(6,0),(3,0)]),
(1.40,[(100,100),(50,95),(25,60),(13,35),(6,10),(3,0)]),
(1.02,[(100,90),(50,70),(25,35),(13,10),(6,10),(3,10)]),
(0.13,[(100,35),(50,20),(25,0),(13,0),(6,0),(3,0)]),
(0.09,[(100,25),(50,25),(25,25),(13,10),(6,0),(3,0)]),
(1.07,[(100,100),(50,45),(25,10),(13,10),(6,0),(3,0)]),
(0.31,[(100,50),(50,20),(25,20),(13,0),(6,0),(3,0)]),
(1.04,[(100,90),(50,70),(25,45),(13,10),(6,0),(3,10)]),
(0.40,[(100,60),(50,20),(25,0),(13,0),(6,0),(3,0)]),
(1.87,[(100,100),(50,100),(25,90),(13,80),(6,60),(3,20)]),
(0.91,[(100,70),(50,70),(25,70),(13,0),(6,35),(3,35)])
]


results2 = [
#(5,[(50,100),(12.5,100),(3.125,100),(0.7813,100),(0.1953,100),(0.0488,100)]),
(3.5,[(50,100),(12.5,70),(3.125,100),(0.7813,100),(0.1953,100),(0.0488,70)]),
#(3,[(50,100),(12.5,100),(3.125,100),(0.7813,100),(0.1953,100),(0.0488,0)]),

]

plt.xlabel('Concentration')