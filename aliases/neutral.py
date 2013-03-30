"""
   This file:     neutral.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Class AliasNeutral wraps aliases and their associated SMILES 
   notations for neutral molecule.

   Copyright (C) 2010  Axel Drefahl

   This file is part of the CurlySMILES package.

   The CurlySMILES package is free software: you can redistribute it
   and/or modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation, either version 3 of
   the License, or (at your option) any later version.

   The CurlySMILES package is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the CurlySMILES package.
   If not, see <http://www.gnu.org/licenses/>.
"""
class AliasNeutral:

   def __init__(self,assignDict=1):      
      self.dictNeu = None
      if assignDict == 1:
         self.initDict()      

   def getSmiles(self,sAlias):
      if self.dictNeu.has_key(sAlias):
         return self.dictNeu[sAlias]
      else:
         return None

   def getDict(self):
      return self.dictNeu

   def initDict(self):

      self.dictNeu = {

         'Ad':         'C1C(C2)CC(CC2C3)CC31', # adamantane

         # 1-naphthylthiourea [86-88-4]
         'ANTU':       'NC(=S)Nc1cccc2ccccc21',

         # 1,1'-(azodicarbonyl)dipiperidine [10465-81-3]
         'ADDP':       'C1CCCCN1C(=O)N=NC(=O)N2CCCCC2',

         # (+/-)-BINOL [602-09-5]
         'BINOL':      'Oc1ccc2ccccc2c1-c3c(O)ccc4ccccc43',

         # 3,3',4,4'-benzophenonetetracarboxylic dianhydride [2421-28-5]
         'BTDA':       'O=C1OC(=O)c2ccc(cc21)C(=O)c3cc4C(=O)OC(=O)c4cc3',
         
         'Bz':         'c1ccccc1',   # benzene
         'CO':         '[C]=O',      # carbon monoxide
         'CO2':        'O=C=O',      # carbon dioxide
         'O3':         'O=[O+][O-]', # ozone

         'TEOS':       'CCO[Si](OCC)(OCC)OCC', # tetraethoxysilane
         'TMOS':       'CO[Si](OC)(OC)OC'      # tetramethoxysilane         

      }

      




