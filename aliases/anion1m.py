"""
   This file:     anion1m.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Class AliasAnion1p wraps aliases and their associated SMILES 
   notations for anions with charge 1-.

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
class AliasAnion1m:

   def __init__(self,assignDict=1):      
      self.dictA1m = None
      if assignDict == 1:
         self.initDict()      

   def getSmiles(self,sAlias):
      if self.dictA1m.has_key(sAlias):
         return self.dictA1m[sAlias]
      else:
         return None

   def getDict(self):
      return self.dictA1m

   def initDict(self):

      self.dictA1m = {
         'Ac(1-)':         '[O-]C(=O)C',

         # acesulfamate (Chem. Commun. 2004, 630-631)
         'Ace(1-)':        'O=S1(=O)[N-]C(=O)C=C(C)O1',


         'AsF6(1-)':       'F[As](F)(F)(F)(F)F{CHe=-1}',
         
         'BF4(1-)':        'F[B](F)(F)[F-]',
         'BPh4(1-)':       'c1ccccc1[B](c2ccccc2)(c3ccccc3)c4ccccc4{CHe=-1}',

         'Br3(1-)':        'Br[Br-]Br', # tribromide


         # cyclopentadienide anion
         'C5H5(1-)':       'c1cccc1{!re=-1}',
         
         # cycloheptatrienyl anion
         'C7H7(1-)':       'C1=CC=CC=C[CH-]1',

         'ClO4(1-)':       'O=[Cl](=O)(=O)[O-]',

         # trifluoromethanesulfonate
         'CF3SO3(1-)':     '[O-]S(=O)(=O)C(F)(F)F',
         
         'CTf3(1-)':\
           'FC(F)(F)S(=O)(=O)[C-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F',        

         # dicyanamide
         'DCA(1-)':        'N#C[N-]C#N',

         # alkylsulfates
         'MeSO4(1-)':      'COS(=O)(=O)[O-]',
         'EtSO4(1-)':      'CCOS(=O)(=O)[O-]',
         'PrSO4(1-)':      'CCCOS(=O)(=O)[O-]',
         'BuSO4(1-)':      'CCCCOS(=O)(=O)[O-]',
         'PeSO4(1-)':      'CCCCCOS(=O)(=O)[O-]',
         'HxSO4(1-)':      'CCCCCCOS(=O)(=O)[O-]',
         'HpSO4(1-)':      'CCCCCCCOS(=O)(=O)[O-]',         
         'OcSO4(1-)':      'CCCCCCCCOS(=O)(=O)[O-]',

         # alkylsulfonates
         'MeSO3(1-)':      'CS(=O)(=O)[O-]',
         'EtSO3(1-)':      'CCS(=O)(=O)[O-]',
         'PrSO3(1-)':      'CCCS(=O)(=O)[O-]',
         'BuSO3(1-)':      'CCCCS(=O)(=O)[O-]',
         'PeSO3(1-)':      'CCCCCS(=O)(=O)[O-]',
         'HxSO3(1-)':      'CCCCCCS(=O)(=O)[O-]',
         'HpSO3(1-)':      'CCCCCCCS(=O)(=O)[O-]',         
         'OcSO3(1-)':      'CCCCCCCCS(=O)(=O)[O-]',
         

         'NO2(1-)':        'O=[N]=O{-1}',  # nitrite
         'NO3(1-)':        '[O-]N(=O)=O)', # nitrate
         'NTf2(1-)':       'FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F',


         'OTf(1-)':     '[O-]S(=O)(=O)C(F)(F)F',

         # tris(pentafluoroethyl)trifluorophosphate
         'PF3(C2F5)3(1-)': \
           '[F-][P](F)(F)(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',

         'PF6(1-)':        'F[P](F)(F)(F)(F)[F-]',
         'AsF6(1-)':       'F[As](F)(F)(F)(F)[F-]',         
         'SbF6(1-)':       'F[Sb](F)(F)(F)(F)[F-]',

         'HSO3(1-)':       'O=[S](O)[O-]',  # hydrogensulfite
         'HSeO3(1-)':      'O=[Se](O)[O-]', # hydrogenselenite

         
         # sacharinate (Chem. Commun. 2004, 630-631)
         'Sac(1-)':        'O=S1(=O)[N-]C(=O)c2ccccc21',

         # thiocyanate
         'SCN(1-)':        'S=C=[N-]',

         # trifluoroacetate
         'TFA(1-)':        'FC(F)(F)C(=O)[O-]',   

         # tosylate
         'Tos(1-)':        '[O-]S(=O)(=O)c1ccc(C)cc1',

         # 2,2,2-trifluoro-N-(trifluoromethylsulfonyl)acetamide
         # (Chem. Commun. 2002, 1726-1727)
         'TSAC(1-)':        'FC(F)(F)S(=O)(=O)[N-]C(=O)C(F)(F)F'         
 
      }

      




