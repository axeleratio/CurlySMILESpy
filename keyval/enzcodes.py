"""
   This file:     enzcodes.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      www.axeleratio.com/axel/axel.htm

   Class EnzymeCodes wraps short codes for enzymes with
   their associated scientific terms.
   To be used as entry in annotation dictionary: {AMenz=...}
   (for example, with operational annotation marker (OPAM):  AM = +Y )

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
class EnzymeCodes:

   def __init__(self,assignDict=1):      
      self.dictEnz = None
      if assignDict == 1:
         self.initDict()      

   def getSciTerm(self,sAlias):
      if self.dictEnz.has_key(sCode):
         return self.dictEnz[sCode]
      else:
         return None

   def getDict(self):
      return self.dictEnz

   def initDict(self):

      self.dictEnz = {

         'SOD':     'superoxide dismutase',
         'TLL':     'Thermomyces Lanuginosa Lipase',
      }
"""
    Notice that this dictionary is growing with future
    CurlySMILES versions as new acronyms/short-notations
    are found in common use in publications regarding
    molecule-enzyme conjugates.
"""




