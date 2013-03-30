"""
   This file:     sec_anion1m.py
   Last modified: October 21, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm   

   Class AliasAnion1m refers secondary aliases for anions with
   charge 1- to their primary aliases in the corresponding module
   in /aliases.

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
      self.dictA1p = None
      if assignDict == 1:
         self.initDict()      

   def getSmiles(self,sAlias):
      if self.dictA1p.has_key(sAlias):
         return self.dictA1p[sAlias]
      else:
         return None

   def getDict(self):
      return self.dictA1p

   def initDict(self):

      self.dictA1p = {

         'Tf2N(1-)':       'NTf2(1-)',
         'TFSI(1-)':       'NTf2(1-)'         

         }

      




