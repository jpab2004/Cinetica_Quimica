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
    'H2': {'name-en': 'Hydrogen gas', 'name-pt': 'Gás hidrogênio', 'radii-empirical': 50, 'radii-calculated': 106, 'mass': 2.016, 'reagents-chance': [[1, 1, .05]], 'color': HEX2VEC('#32b82e')},
    'OH': {'name-en': 'Hydroxy group', 'name-pt': 'Hidroxila', 'radii-empirical': 97, 'radii-calculated': 101, 'mass': 17.00684, 'reagents-chance': [[1, 8, .25]], 'color': HEX2VEC('#ff9999')},
    'H2O': {'name-en': 'Water', 'name-pt': 'Água', 'radii-empirical': 75, 'radii-calculated': 150, 'mass': 18.01528, 'reagents-chance': [['H2', 8, .01], ['OH', 1, .1]], 'color': HEX2VEC('#0000ff')}
}



# Saving file
fileName = 'Assets/moleculesHashTable.pickle'

with open(fileName, 'wb') as f:
    dump(molecules, f, protocol=HIGHEST_PROTOCOL)

with open(fileName, 'rb') as f:
    molecules2 = load(f)

print(molecules == molecules2)