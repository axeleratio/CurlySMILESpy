"""
   This file:     csm_curlyann.py
   Last modified: November 9, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Python module csm_curlyann implements a class for managing
   an annotation enclosed in curly braces.

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
import csm_dataface
class CurlyAnnotation:
    
   def __init__(self,oDataFace=None,sContent=None):

      self.oDataFace = oDataFace
      self.sContent  = sContent # content inside curly braces
      self.sAM     = None # annotation marker (one- or two-char)
      self.sAMtype = None # annotation marker type   
      self.dictAnn = None # annotation dictionary      
      self.lstErrors = [] # list of strings with reported error 

   #===================================================================#
   # PARSE annotation                                                         #
   #===================================================================#

   """------------------------------------------------------------------
      parse: parse self.sContent and assign self.sAM and self.dictAnn
      return: self.lstErrors
   """
   def parse(self):

      self.dictAnn = {}

      nLenCont = len(self.sContent)
      if nLenCont < 3:
         self.sAM = self.sContent        
      else:
         self.sAM   = self.sContent[0:2]  
         sDict = self.sContent[2:]
         if len(sDict) > 2:
            self.parse_dict(sDict)
         else:
            sMsg = "parse: incomplete annotation dictionary: '%s'" %\
                   sDict
            self.lstErrors.append(sMsg)
   
      self.sAMtype = self.oDataFace.marker_type(self.sAM)
      if self.sAMtype == None:
         sMsg = "parse: unknown descriptor or marker: '%s'" % self.sAM
         self.lstErrors.append(sMsg)
   
      return self.lstErrors


   """------------------------------------------------------------------
      parse: parse annotation dictionary sDict and map its content
             into self.dictAnn
   """
   def parse_dict(self,sDict):

      # get list with dictionary entries, each a pair of a key and
      # a value separated by an equal sign 
      lstParts = self.split_at_top_level(sDict,';')

      for sPart in lstParts:
         pair = sPart.split('=',1) # split at first equal sign
         if len(pair) != 2:
            sMsg  = "parse_dict: invalid format for annotation "
            sMsg += "entry: '%s'" % sPart
            self.lstErrors.append(sMsg)
         else:
            sKey = pair[0]
            sVal = pair[1]
            if len(sKey) < 1:
               sMsg  = "parse_dict: key of entry '%s' is empty" % sPart
               self.lstErrors.append(sMsg)
               return
            if len(sVal) < 1:
               sMsg  = "parse_dict: value of entry '%s' is empty" % sPart
               self.lstErrors.append(sMsg)
               return
            if sKey[0] != '$' and not self.oDataFace.is_anndict_key(sKey):
               sMsg  = "parse_dict: key '%s' not defined" % sKey
               self.lstErrors.append(sMsg)
               return               
            self.dictAnn[sKey] = sVal

   """
      split_at_top_level: split given string at top level by cSplit,
                          i.e do not split if cSplit is curly-nested  
      return: lstParts
   """
   def split_at_top_level(self,sStr,cSplit):
      lstParts        = []
      sSub          = ''
      iCurlyDepth = 0 # curly brackets nesting depth
      for cChar in sStr:
         if cChar == cSplit:
            if iCurlyDepth == 0: # split-point reached
               lstParts.append(sSub)
               sSub = ''
               continue
         sSub += cChar         
         if cChar == '{':
            iCurlyDepth += 1
         elif cChar == '}': 
            iCurlyDepth -= 1

      if iCurlyDepth != 0:
         sErr = 'csm_curlycode.splitAtTopLevel [%s]: ' % self.sIndexes
         sErr = 'unbalanced curly entry: "%s"' % sStr
         self.lstErrors.append(sErr)      
      else:
         lstParts.append(sSub) # append final part

      return lstParts

   #===================================================================#
   # VERIFY keys in annotation dictionary                              #
   #===================================================================#
   """
      split_at_top_level: split given string at top level by cSplit,
                          i.e do not split if cSplit is curly-nested  
      return: (lstClient,lstUnknown)
               lstClient  = list of client keys (starting with $)
               lstUnknown = list of unknown keys
   """
   def verify_keys(self):
      lstClient  = []
      lstUnknown = []

      for sKey in self.dictAnn.keys():
         if sKey[0] == '$':
            lstClient.append(sKey)
         elif not self.oDataFace.is_anndict_key(sKey):
            lstUnknown.append(sKey)

      return (lstClient,lstUnknown)  

   #===================================================================#
   # ACCESS member variables from application                          #
   #===================================================================#
   def content(self):                 return self.sContent
   def annotation_marker(self):       return self.sAM
   def annotation_marker_type(self):  return self.sAMtype   
   def annotation_dict(self):         return self.dictAnn 
   def msgs_err(self):                return self.lstErrors 

   """------------------------------------------------------------------
      entry: get content as marker/dict pair  
      return: [self.sAM, self.dictAnn]
   """
   def entry(self):
      return [self.sAM, self.dictAnn]
   
   #===================================================================#
   # QUERIES                                                           #
   #===================================================================#
   def is_stereodescr(self):
      return self.oDataFace.is_stereodescriptor(self.sAM) 
   def is_open_bond(self):
      return self.oDataFace.is_bond_symbol(self.sAM)
   def is_GEA(self): #group environment annotation
      return self.oDataFace.is_GEAM(self.sAM)
   def is_MDA(self): # molecular detail annotation?
      return self.oDataFace.is_MDAM(self.sAM)
   def is_OPA(self): # operational annotation?
      return self.oDataFace.is_OPAM(self.sAM)
   def is_SSA(self): # state and shape annotation?
      return self.oDataFace.is_SSAM(self.sAM)      
   def is_MIA(self): # miscellaneous interest annotation?
      return self.oDataFace.is_MIAM(self.sAM)
   
