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
molecules = {
    'H2': {'name-en': 'Hydrogen gas', 'name-pt': 'Gás hidrogênio', 'radii-empirical': 50, 'radii-calculated': 106, 'mass': 2.016, 'atom-types': [1, 1], 'color': HEX2VEC('#f55f73')}
}



# Saving file
fileName = 'Assets/moleculesHashTable.pickle'

with open(fileName, 'wb') as f:
    dump(molecules, f, protocol=HIGHEST_PROTOCOL)

with open(fileName, 'rb') as f:
    molecules2 = load(f)

print(molecules == molecules2)