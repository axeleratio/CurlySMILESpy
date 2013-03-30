"""
   This file:     procodes.py
   Last modified: October 21, 2010
   Package:       CurlySMILES Version 1.0.0
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com

   Class ProteinCodes wraps short codes for proteins with
   their associated scientific terms
   To be used as entry in annotation dictionary: {AMpro=...}
   
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
class ProteinCodes:

   def __init__(self,assignDict=1):      
      self.dictPro = None
      if assignDict == 1:
         self.initDict()      

   def getSciTerm(self,sAlias):
      if self.dictPro.has_key(sCode):
         return self.dictPro[sCode]
      else:
         return None

   def getDict(self):
      return self.dictPro

   def initDict(self):

      self.dictPro = {

         'BSA':     'bovine serum albumin',
      }
"""
    Notice that this dictionary is growing with future
    CurlySMILES versions as new acronyms/short-notations
    are found in common use in publications regarding
    molecule-protein conjugates.
"""
      




