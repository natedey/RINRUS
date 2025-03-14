import os, sys, re
from numpy import *
from read_write_pdb import *
from copy import *
import argparse
import operator
import pandas as pd
import numpy as np


Atoms = {'H' : [1,  1.00794],
         'HE': [2,  4.002602],
         'LI': [3,  6.941],
         'BE': [4,  9.012182],
         'B' : [5,  10.811],
         'C' : [6,  12.0107],
         'N' : [7,  14.00674],
         'O' : [8,  15.9994], 
         'F' : [9,  18.9984],
         'NE': [10, 20.1797],
         'NA': [11, 22.9897],
         'MG': [12, 24.305],
         'AL': [13, 26.9815],
         'SI': [14, 28.0855],
         'P' : [15, 30.9738],
         'S':  [16, 32.065],
         'CL': [17, 35.453],
         'Ar': [18, 39.962383], 
         'K' : [19, 39.0983],
         'CA': [20, 40.078],
         'SC': [21, 44.955912],
         'TI': [22, 47.867000],
         'V' : [23, 50.9415],
         'CR': [24, 51.9961],
         'MN': [25, 54.938],
         'FE': [26, 55.845],
         'NI': [27, 58.6934],
         'CO': [28, 58.9332],
         'CU': [29, 62.939598],
         'ZN': [30, 65.39],
         'GA': [31, 69.723],
         'GE': [32, 72.64],
         'AS': [33, 74.9216],
         'SE': [34, 78.96],
         'BR': [35, 79.904],
         'KR': [36, 83.8],
         'RB': [37, 85.4678],
         'SR': [38, 87.62],
         'Y' : [39, 88.9059],
         'ZR': [40, 91.224],
         'NB': [41, 92.9064],
         'MO': [42, 95.94],
         'TC': [43, 98],
         'RU': [44, 101.07],
         'RH': [45, 102.9055],
         'PD': [46, 106.42],
         'AG': [47, 107.8682],
         'CD': [48, 112.411],
         'IN': [49, 114.818],
         'SN': [50, 118.71],
         'Sb': [51, 121.76],        #Antimony            Sb  51
         'I' : [53, 126.9045],      #Iodine          I   53
         'Te': [52, 127.6],         #Tellurium           Te  52
         'Xe': [54, 131.293],       #Xenon               Xe  54
         'Cs': [55, 132.9055],      #Cesium          Cs  55
         'Ba': [56, 137.327],       #Barium              Ba  56
         'La': [57, 138.9055],      #Lanthanum       La  57
         'Ce': [58, 140.116],       #Cerium              Ce  58
         'Pr': [59, 140.9077],      #Praseodymium    Pr  59
         'Nd': [60, 144.24],        #Neodymium           Nd  60
         'Pm': [61, 145],           #Promethium              Pm  61
         'Sm': [62, 150.36],        #Samarium            Sm  62
         'Eu': [63, 151.964],       #Europium            Eu  63
         'Gd': [64, 157.25],        #Gadolinium          Gd  64
         'Tb': [65, 158.9253],      #Terbium         Tb  65
         'Dy': [66, 162.5],         #Dysprosium          Dy  66
         'Ho': [67, 164.9303],      #Holmium         Ho  67
         'Er': [68, 167.259],       #Erbium              Er  68
         'Tm': [69, 168.9342],      #Thulium         Tm  69
         'Yb': [70, 173.04],        #Ytterbium           Yb  70
         'Lu': [71, 174.967],       #Lutetium            Lu  71
         'Hf': [72, 178.49],        #Hafnium             Hf  72
         'Ta': [73, 180.9479],      #Tantalum        Ta  73
         'W' : [74, 183.84],        #Tungsten            W   74
         'Re': [75, 186.207],       #Rhenium             Re  75
         'Os': [76, 190.23],        #Osmium              Os  76
         'Ir': [77, 192.217],       #Iridium             Ir  77
         'Pt': [78, 195.078],       #Platinum            Pt  78
         'Au': [79, 196.9665],      #Gold            Au  79
         'Hg': [80, 200.59],        #Mercury             Hg  80
         'Tl': [81, 204.3833],      #Thallium        Tl  81
         'Pb': [82, 207.2],         #Lead                Pb  82
         'Bi': [83, 208.9804],      #Bismuth         Bi  83
         'Po': [84, 209],           #Polonium                Po  84
         'At': [85, 210],           #Astatine                At  85
         'Rn': [86, 222],           #Radon                   Rn  86
         'Fr': [87, 223],           #Francium                Fr  87
         'Ra': [88, 226],           #Radium                  Ra  88
         'Ac': [89, 227],           #Actinium                Ac  89
         'Pa': [91, 231.0359],      # Protactinium    Pa  91
         'Th': [90, 232.0381],      # Thorium         Th  90
         'Np': [93, 237],           # Neptunium               Np  93
         'U' : [92, 238.0289],      # Uranium         U   92
         'Am': [95, 243],           # Americium               Am  95
         'Pu': [94, 244],           # Plutonium               Pu  94
         'Cm': [96, 247],           # Curium                  Cm  96
         'Bk': [97, 247],           # Berkelium               Bk  97
         'Cf': [98, 251],           # Californium             Cf  98
         'Es': [99, 252],           # Einsteinium             Es  99
         'Fm': [100,    257],       #Fermium                 Fm  100 
         'Md': [101,    258],       #Mendelevium             Md  101 
         'No': [102,    259],       #Nobelium                No  102 
         'Rf': [104,    261],       #Rutherfordium           Rf  104 
         'Lr': [103,    262],       #Lawrencium              Lr  103 
         'Db': [105,    262],       #Dubnium                 Db  105 
         'Bh': [107,    264],       #Bohrium                 Bh  107 
         'Sg': [106,    266],       #Seaborgium              Sg  106 
         'Mt': [109,    268],       #Meitnerium              Mt  109 
         'Rg': [111,    272],       #Roentgenium             Rg  111 
         'Hs': [108,    277],       #Hassium                 Hs  108 
         'Ds': [110,    0],         #Darmstadtium                Ds  110 
         'Cn': [112,    0],         #Copernicium                 Cn  112 
         'Nh': [113,    0],         #Nihonium                    Nh  113 
         'Fl': [114,    0],         #Flerovium                   Fl  114 
         'Mc': [115,    0],         #Moscovium                   Mc  115 
         'Lv': [116,    0],         #Livermorium                 Lv  116 
         'Ts': [117,    0],         #Tennessine                  Ts  117 
         'Og': [118,    0]          #Oganesson                   Og  118 
         }


def center_of_mass(mass,coord):
    mass = mass.astype(float)[:,newaxis]
    coord = coord.astype(float)
    return sum(coord*mass,axis=0)/sum(mass)


def avg_coord(coord):      ## average coordinates of a residue
    return mean(coord,axis=0)


def calc_dist(a,b):  ## a and b are xyz coordinates array
    dist = sqrt( (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 )
    return float(dist)


