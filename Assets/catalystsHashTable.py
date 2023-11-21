from pickle import dump, load, HIGHEST_PROTOCOL
from vpython import vector



# Functions
def RGB2VEC(r, g, b): 
    return vector(r/255, g/255, b/255)

def HEX2VEC(hex):
    hex = hex.lstrip('#')
    r, g, b = tuple(int(hex[i:i+2], 16) for i in (0, 2, 4))
    return RGB2VEC(r, g, b)



# Dicionaries
catalysts = {
    'Fe': {1: [[300, .01], [300, 1.5]]}
}



# Saving file
fileName = 'Assets/catalystsHashTable.pickle'

with open(fileName, 'wb') as f:
    dump(catalysts, f, protocol=HIGHEST_PROTOCOL)

with open(fileName, 'rb') as f:
    molecules2 = load(f)

print(catalysts == molecules2)