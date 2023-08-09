from pickle import dump, load
from vpython import vector

def RGB2VEC(r, g, b): 
    return vector(r/255, g/255, b/255)

def HEX2VEC(hex):
    hex = hex.lstrip('#')
    r, g, b = tuple(int(hex[i:i+2], 16) for i in (0, 2, 4))
    return RGB2VEC(r, g, b)

elements = {
    1: { 'name': 'Hydrogen', 'symbol': 'H', 'radius': 25, 'mass': 1.00797, 'color': HEX2VEC('#ffffff') },
    2: { 'name': 'Helium', 'symbol': 'He', 'radius': 120, 'mass': 4.00260, 'color': HEX2VEC('#d9ffff') },
    3: { 'name': 'Lithium', 'symbol': 'Li', 'radius': 145, 'mass': 6.941, 'color': HEX2VEC('#cc80ff') },
    4: { 'name': 'Beryllium', 'symbol': 'Be', 'radius': 105, 'mass': 9.01218, 'color': HEX2VEC('#c2ff00') },
    5: { 'name': 'Boron', 'symbol': 'B', 'radius': 85, 'mass': 10.81, 'color': HEX2VEC('#ffb5b5') },
    6: { 'name': 'Carbon', 'symbol': 'C', 'radius': 70, 'mass': 12.011, 'color': HEX2VEC('#909090') },
    7: { 'name': 'Nitrogen', 'symbol': 'N', 'radius': 65, 'mass': 14.0067, 'color': HEX2VEC('#3050f8') },
    8: { 'name': 'Oxygen', 'symbol': 'O', 'radius': 60, 'mass': 15.9994, 'color': HEX2VEC('#ff0d0d') },
    9: { 'name': 'Fluorine', 'symbol': 'F', 'radius': 50, 'mass': 18.998403, 'color': HEX2VEC('#90e050') },
    10: { 'name': 'Neon', 'symbol': 'Ne', 'radius': 160, 'mass': 20.179, 'color': HEX2VEC('#b3e3f5') },
    11: { 'name': 'Sodium', 'symbol': 'Na', 'radius': 180, 'mass': 22.98977, 'color': HEX2VEC('#ab5cf2') }
}

for e in elements:
    print(elements[e]['color'])