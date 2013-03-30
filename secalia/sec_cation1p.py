"""
   This file:     sec_cation1p.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Class AliasCation1m refers secondary aliases for cations with
   charge 1+ to their primary aliases in the corresponding module
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
class AliasCation1p:

   def __init__(self,assignDict=1):      
      self.dictC1p = None
      if assignDict == 1:
         self.initDict()      

   def getSmiles(self,sAlias):
      if self.dictC1p.has_key(sAlias):
         return self.dictC1p[sAlias]
      else:
         return None

   def getDict(self):
      return self.dictC1p

   def initDict(self):

      self.dictC1p = {

         # imidazolium
         'emim(1+)':       'C2mim(1+)',
         'bmim(1+)':       'C4mim(1+)',
         'omim(1+)':       'C8mim(1+)',
         }

      




