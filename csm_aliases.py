"""
   This file:     csm_aliases.py
   Last modified: October 21, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Python module csm_aliases manages primary and secondary
   aliases, which are replaced by a component notation (SFN,
   Composite, SMILES or annotated SMILES code), when a user
   notation is turned into a work notation.
   The central method is compnt_notation(self, sAlias)
   to get the component notation for an alias.

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
import sys, os

class AliasNotations:
    
   def __init__(self,sCsmpyDir,assignDict=1):
      """
         dictionaries with aliases:
            dictPrimAliases with primary aliases
            dictSecAliases with secondary aliases (for a secondary alias
                                 a corresponding primary one must exist)
            dictClientAliases with client-provided aliases                   
      """
      self.sCsmpyDir = sCsmpyDir
      self.dictPrimAliases = None
      self.dictSecAliases = None
      if assignDict == 1:
         self.initDict()

   def initDict(self):
       self.dictPrimAliases = self.load_prim_aliases()
       self.dictSecAliases  = self.load_sec_aliases()

   #===================================================================#
   # LOAD aliases                                                      #
   #===================================================================#
   """------------------------------------------------------------------   
      loadPrimAliases:
         load primary aliases from dictionaries in python modules 
         located in directory aliases (relative to directory of this
         module, which is expected to be in subdirectory csm/py);
         
      return: dictAliases = {sAliasGroup: dictGroup,...}
                 sAliasGroup = alias group name such as neutral,
                                     cation1p, anion1p, etc., equal to
                                     module name in directory aliases;
                 dictGroup = {sPrimAlias:sSmiles,...}              
   """
   def load_prim_aliases(self):
            
      (sDirAliases,lstPyMod) = self.get_module_paths('aliases')
      
      # get all dictionaries with primary aliases
      dictPrimAliases = {}
      code0 = "sys.path.append('%s')" % sDirAliases
      exec code0
      for sPyMod in lstPyMod:
         lstParts = sPyMod.split(os.sep) 
         sAliasGroup = lstParts[-1][0:-3]
         code1 = "import %s" % sAliasGroup 
         exec code1
         sClassName = sAliasGroup[0].upper() + sAliasGroup[1:]
         sClassName = 'Alias' + sClassName
         code2 = "oDict = %s.%s()" % (sAliasGroup,sClassName)
         exec code2         
         dictPrimAliases[sAliasGroup] = oDict.getDict()
         del oDict        

      return dictPrimAliases

   """------------------------------------------------------------------
      load_sec_aliases:
         load secondary aliases from dictionaries in python modules 
         located in directory secalia (relative to directory of this 
         this module, which is expected to be in subdirectory csm/py);
         
      return: dictAliases = {sAliasGroup: dictGroup,...}
                 sAliasGroup = alias group name such as neutral,
                                     cation1p, anion1p, etc., equal to
                                     module name in directory aliases;
                 dictGroup = {sSecAlias:sPrimAlias,...}              
   """
   def load_sec_aliases(self):
            
      (sDirAliases,lstPyMod) = self.get_module_paths('secalia')
      
      # get all dictionaries with secondary aliases
      dictSecAliases = {}
      code0 = "sys.path.append('%s')" % sDirAliases
      exec code0
      for sPyMod in lstPyMod:
         lstParts = sPyMod.split(os.sep) 
         sAliasGroup = lstParts[-1][0:-3]

         # take only modules having a name starting with 'sec_'
         if cmp(sAliasGroup[0:4],'sec_') == 0:
            sAliasGroup = sAliasGroup[4:]
         else:
            continue
         
         code1 = "import sec_%s" % sAliasGroup 
         exec code1
         sClassName = sAliasGroup[0].upper() + sAliasGroup[1:]
         sClassName = 'Alias' + sClassName
         code2 = "oDict = sec_%s.%s()" % (sAliasGroup,sClassName)
         exec code2         
         dictSecAliases[sAliasGroup] = oDict.getDict()
         del oDict        

      return dictSecAliases

   """------------------------------------------------------------------
      get_module_paths: find and list absolute path for
                      each alias module either in sub-directory
      sSubdir = 'aliases' or 'secalia'               

      return: (sDirAliases,lstPyMod), where
               sDirAliases is absolute path to subdirectory, and
               lstPyMod is list of absolute paths to module files
   """
   def get_module_paths(self,sSubdir):

      # absolute path to aliases directory
      sDirAliases =  self.sCsmpyDir + os.sep + sSubdir 

      # get names of aliases modules
      lstPyMod = []
      lstFiles = []
      lstFiles = os.listdir(sDirAliases)
      for sFile in lstFiles:
         if len(sFile) < 6 or sFile[-3:] != '.py':
            continue
         sCompletePath = sDirAliases + os.sep + sFile
         if os.path.isfile(sCompletePath):
            lstPyMod.append(sCompletePath)

      return (sDirAliases,lstPyMod)
  
   #===================================================================#
   # LOOK-UP alias (and group-id)                                      #
   #===================================================================#
   """------------------------------------------------------------------
      compnt_notation: look up alias        
      return: string with component notation; or 
              None, if alias not found
   """
   def compnt_notation(self, sAlias):
      (sCompntNotation,sGroupID) = \
         self.compnt_notation_and_groupid(sAlias)
      return sCompntNotation 
  
   """------------------------------------------------------------------
      compnt_notation_and_groupid: look up alias 

      return:  (sCompntNotation,sGroupId) or (None,None), if not found
                sCompntNotation  = notation for alias replacement 
                sGroupId = alias group name such as neutral, cation1p,
                           anion1p,etc.
   """
   def compnt_notation_and_groupid(self,sAlias):
      sCompntNotation = None
      sGroupId = None

      # Primary alias first ...
      (sCompntNotation,sGroupId) = self.lookup_as_prim_alias(sAlias)
      # ... if not found ...
      if sCompntNotation == None: # look up as secondary alias
         (sPrimAlias,sCompntNotation,sGroupId) = \
            self.lookup_as_sec_alias(sAlias)                           

      return (sCompntNotation,sGroupId)

   """------------------------------------------------------------------
      lookup_as_prim_alias:
      return: (sCompntNotation,sGroupId) for primary alias or
              (None,None) if not found
   """      
   def lookup_as_prim_alias(self,sPrimAlias):
      for sGroupId in self.dictPrimAliases:
         dict = self.dictPrimAliases[sGroupId]
         if dict.has_key(sPrimAlias):
            return (dict[sPrimAlias],sGroupId)
      return (None,None)

   """------------------------------------------------------------------
      lookup_as_prim_alias_by_groupid:
      return: sCompntNotation for primary alias or None if not found
   """
   def lookup_as_prim_alias_by_groupid(self,sPrimAlias,sGroupId):
      if self.dictPrimAliases.has_key(sGroupId):
         dict = self.dictPrimAliases[sGroupId]
         if dict.has_key(sPrimAlias):
            return dict[sPrimAlias]
      else:   
         return None

   """------------------------------------------------------------------
      lookup_as_sec_alias:
      return: (sPrimAlias, sCompntNotation,sGroupId) for secondary alias 
              or (None,None,None) if not found
   """      
   def lookup_as_sec_alias(self,sSecAlias):
      for sGroupId in self.dictSecAliases:
         dict = self.dictSecAliases[sGroupId]
         if dict.has_key(sSecAlias):
            sPrimAlias = dict[sSecAlias]
            sCompntNotation = \
               self.lookup_as_prim_alias_by_groupid(sPrimAlias,sGroupId)
            if sCompntNotation != None:
               return (sPrimAlias,sCompntNotation,sGroupId)
      return (None,None,None)

   #===================================================================#
   # MAKE alias dictionary containing primary and secondary aliases    #
   #===================================================================#
   """------------------------------------------------------------------
       makeAliasDict: check consistency of alias-alias and alias-groupid
                      relations and make dictionary that has both 
                      primary and secondary aliases as key, while value
                      is the corresponding primary alias (if key is a
                      primary alias then key and value are the same)

                      NOTE: this method is for use during development
                            and extension of alias dictionaries

       return:  (lstAmbig,dictAliases)
                lstAmbig = list of lines, each line reporting an
                           ambiguity (multiply used alias name)
                           empty if no ambiguities
                dictAliases: {sAlias: sPrimAlias,...}
                              sAlias = primary or secondary alias
                              sPrimAlias = primary alias corresponding
                                           to sAlias
       Note: client aliases are not considered here         
   """
   def makeAliasDict(self):
      lstAmbig = []
      dictAliases = {}

      # primary aliases
      dictPrimGroupId = {} # dict with first encountered group id
      for sPrimGroupId in self.dictPrimAliases:
         dictPrim = self.dictPrimAliases[sPrimGroupId]
         for sPrimAlias in dictPrim:
            if dictAliases.has_key(sPrimAlias):
               sLine = '"%s" with two group ids: "%s" and "%s"' % \
                 (sPrimAlias,sPrimGroupId, dictPrimGroupId[sPrimAlias])
               sLine += ' (both for primary alias)'
               lstAmbig.append(sLine)
            else:
               dictPrimGroupId[sPrimAlias] = sPrimGroupId 
               dictAliases[sPrimAlias] = sPrimAlias

      # secondary aliases   
      dictSecGroupId = {} # dict with first encountered group id
      for sSecGroupId in self.dictSecAliases:
         dictSec = self.dictSecAliases[sSecGroupId]
         for sSecAlias in dictSec:
            sPrimAliasCorresp = dictSec[sSecAlias]

            # first, check if sec. alias was already used as prim. alias       
            if dictAliases.has_key(sSecAlias):
               sLine  = 'sec. alias "%s" ' % sSecAlias
               sLine += 'with group id "%s" conflicts ' % sSecGroupId
               sLine += 'with same-name prim. alias of group "%s"' % \
                     dictPrimGroupId[sSecAlias]  
               lstAmbig.append(sLine)
               continue
            # also make sure the corresp. prim. alias exists
            elif not dictAliases.has_key(sPrimAliasCorresp):
               sLine  = 'sec. alias "%s" ' % sSecAlias
               sLine += 'with group id "%s" ' % sSecGroupId               
               sLine += 'has no corresponding prim. alias '
               sLine += 'named "%s"' % sPrimAliasCorresp
               lstAmbig.append(sLine)
               continue
            else:
               # also make sure prim. and sec. share same group id 
               (sSmiles,sGroupIdCorresp) = \
                  self.lookupAsPrimAlias(sPrimAliasCorresp)
               if cmp(sSecGroupId,sGroupIdCorresp) != 0:
                  sLine  = 'group id mismatch for sec. alias '
                  sLine += '"%s" in group "%s": ' % \
                           (sSecAlias,sSecGroupId)  
                  sLine += 'corresp. prim. alias "%s" ' % \
                           sPrimAliasCorresp
                  sLine += 'is in group "%s"' % sGroupIdCorresp
                  lstAmbig.append(sLine)
                  continue

            # check if sec. alias is used twice   
            if dictSecGroupId.has_key(sSecAlias):
               sLine = '"%s" with two group ids: "%s" and "%s"' % \
                 (sSecAlias,sSecGroupId, dictSecGroupId[sPrimAlias])
               sLine += ' (both for secondary alias)'
               lstAmbig.append(sLine)
            else:
               dictSecGroupId[sSecAlias] = sSecGroupId 
               dictAliases[sSecAlias] = sPrimAliasCorresp
            
      return (lstAmbig,dictAliases)
