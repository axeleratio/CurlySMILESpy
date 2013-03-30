"""
   This file:     anion2m.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Class AliasAnion2p wraps aliases and their associated SMILES 
   notations for anions with charge 2-. 

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
class AliasAnion2m:

   def __init__(self,assignDict=1):      
      self.dictA2m = None
      if assignDict == 1:
         self.initDict()      

   def getSmiles(self,sAlias):
      if self.dictA2m.has_key(sAlias):
         return self.dictA2m[sAlias]
      else:
         return None

   def getDict(self):
      return self.dictA2m

   def initDict(self):

      self.dictA2m = {
         'CO3(2-)':         '[O-]C(=O)[O-]',

          # cyclooctatetraenide dianion
         'C8H8(2-)':        'c1ccccccc1{!re=-2}',
         
         'HPO4(2-)':        '[O-]P(=O)(O)[O-]',
         'S(2-)':           '[S-2]',                      # sulfide         
         'SO3(2-)':         '[O]~[S](~[O])~[O]{CHe=-2}',  # sulfite
         'SO4(2-)':         '[O-]S(=O)(=O)[O-]',          # sulfate
         'SeO3(2-)':        '[O]~[Se](~[O])~[O]{CHe=-2}', # selenite
      }

      




