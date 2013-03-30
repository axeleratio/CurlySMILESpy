"""
   This file:     csm_notation.py
   Last modified: November 8, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Python module csm_notation implements a class for managing
   a CurlySMILES notation.

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
import csm_dataface, csm_annsmi, csm_molform, csm_sfn, csm_curlyann
import re
class Notation:
    
   def __init__(self,oDataFace=None,sUserNotation=None):

      # native and client data
      self.oDataFace = oDataFace # reference to object with native data
      self.dictClientAnn  = {}   # dictionary with client annotations
      self.dictClientAli  = {}   # dictionary with client aliases

      # notations
      self.sUserNotation = sUserNotation
      self.sWorkNotation = None

      # components of work notation
      self.nCompnt    = None # number of components in work notation
      self.lstCompnt  = [] # list of component notation
      self.lstTypes   = [] # list of component types: 'cps', 'sfn', 
                             #    'smi', or 'ali' 
      self.lstObjSmi  = [] # list of SMILES notation objects
      self.lstObjSfn  = [] # list of SFN notation objects
      self.lstMfHill = [] # list with molecular formulae of compnents
                          # (if 'smi' or 'sfn', else None) in Hill format

      # parts of composite notation: sComposite = self.lstCompnt[0]
      self.nCpsParts    = [] # number of composite parts      
      self.lstCpsParts  = [] # list with notations of composite parts
      self.lstCpsTypes  = [] # list of classifiers for composite parts:
                                # 'smi' or 'sfn'
      self.lstObjSmiCps = [] # list of SMILES obj for composite parts
      self.lstObjSfnCps = [] # list of SFN obj for composite parts

      # component-anchored annotations
      self.lstCompntAnn = [] # list of dictionaries with
                               # component-anchored dictionaries
       
      # status information and error reports derived while parsing
      self.lstAccFail = [] # list of strings with access failure notes
      self.lstErrors  = [] # list of strings with reported error 

   #===================================================================#
   # SET data needed in parsing a CurlySMILES notation                 #
   #===================================================================#
   """------------------------------------------------------------------
      set_user_notation: assign or reassign self.sUserNotation 
   """
   def set_user_notation(self,sUserNotation):
      self.__init__(self.oDataFace)
      self.sUserNotation = sUserNotation
      
   """------------------------------------------------------------------
      set_client_annotations: assign self.dictClientAnn 
   """
   def set_client_ann(self,dictClientAnn):
      self.dictClientAnn = dictClientAnn

   """------------------------------------------------------------------
      set_client_aliases: assign self.dictClientAliases 
   """
   def set_client_ali(self,dictClientAli):
      self.dictClientAli = dictClientAli


   #===================================================================#
   # PARSE CurlySMILES notation                                        #
   #===================================================================#

   """------------------------------------------------------------------
      parse: parse CurlySMILES user notation by
             preprocessing: turn user into work notation
             evaluating:

      argument: sUserNotation, to be assigned to self.sUserNotation
                               unless done with self.set_user_notation

      return: error_msgs
   """
   def parse(self, sUserNotation=None):

      # assign user notation unless already done
      if self.sUserNotation == None:
         self.sUserNotation = sUserNotation

      # check for fatal errors of highest level
      if self.sUserNotation == None: # if still None
         sMsg = 'parse: sUserNotation equals None'
         self.lstErrors.append(sMsg)
         return self.lstErrors
      elif len(self.sUserNotation) < 1:         
         return self.lstErrors      
      else:   
         for ch in self.sUserNotation:
            if ch in ['<','>','"',"'"]:
               sMsg  = "Notation.parse: invalid char in notation: '%s'" % ch
               self.lstErrors.append(sMsg)
               return self.lstErrors

      # Preprocess
      sModified = self.replace_client_ann(self.sUserNotation)

      # Eliminate white space
      sModified = self.eliminate_backslashes_and_white_space(sModified)

      # Replace aliases based on dict given by client
      sModified = self.replace_client_ali(sModified)

      # Replace aliases based on dict in /aliases and /secalia
      sModified = self.replace_builtin_ali(sModified)

      # Make work notation and associated parameters
      self.make_work_notation(sModified)

      # loop over components to grap error messages
      iCompnt = 0
      for sType in self.lstTypes:

         lstErr = []
         if cmp(sType,'smi') == 0:
            oSmi = self.lstObjSmi[iCompnt]
            lstErr = oSmi.msgs_err()
            for sErr in lstErr:
               sMsg  = 'Notation.parse: %d. component '\
                       % (iCompnt+1)
               sMsg += '< %s' % sErr
               self.lstErrors.append(sMsg)
               
         elif cmp(sType,'sfn') == 0:            
            oSfn = self.lstObjSfn[iCompnt]
            lstErr = oSfn.msgs_err()

         iCompnt += 1
      
      return self.lstErrors      

   #===================================================================#
   # REPLACE client notations                                          #
   #===================================================================#
   """------------------------------------------------------------------
      replace_client_annotations: modify string with current notation
                                  by replacing client annotations if
      replacement code is available from self.dictClientAnn and/or 
      from inline code (via \\ format).

      return: string with modified notation      
   """
   def replace_client_ann(self,sCurrNotation):

      # get parts separated by backslash pairs;
      # first part is main notation, following parts are dictionary
      # entries (annotation-identifier/annotation-content pairs)
      # to be mapped into dictInlineAnn
      lstParts = sCurrNotation.split('\\\\')
      dictInlineAnn = {}
      
      # build dictInlineAnn, if entries found inline 
      if len(lstParts) > 1:  
         for sPart in lstParts[1:]:
            sPartStripped = sPart.strip()
            pair = sPartStripped.split()
            if len(pair) != 2:
              sMsg  = "replace_client_ann: 2 instead of %d " % len(pair)
              sMsg += "white-space-separated strings expected in '%s'" \
                      % sPartStripped
              self.lstErrors.append(sMsg)
              return lstParts[0]           
            elif len(pair[0]) < 3:
              sMsg  = "replace_client_ann: annotation identifier "
              sMsg += "'%s' too short" % pair[0]
              self.lstErrors.append(sMsg)
              return lstParts[0]                          
         
            dictInlineAnn[pair[0]] = pair[1]

      # replace atom-anchored client-annotations using dictInlineAnn
      # and self.dictClientAnn
      sModified = self.replace_by_marker(lstParts[0],dictInlineAnn,'$*')

      return sModified

   """------------------------------------------------------------------
      replace_by_marker: replace
         (a) atom-anchored client annotation, if sMarker='$*', or
         (b) component-anchored client annotation, if sMarker='$$'
          by finding replacement code in self.dictClientAnn and/or
          dictInlineAnn.

      return: string with modified notation      
   """
   def replace_by_marker(self,sCurrNotation,dictInlineAnn,sMarker):

      # break at each point where client annotation begins,
      # which is at '{$*' or '{$$' 
      pat =re.compile('\{\$(\$|\*)') # '*' or '$' is referred back
      lstParts = pat.split(sCurrNotation)

      # if client annotation is identifier, try to replace
      sModified = lstParts[0]
      j = 0
      cMarker = None  # will hold either '*' or '$'
      for sPart in lstParts[1:]:
         j += 1
         if j%2 != 0: # when j odd sPart is '*' or '$'
            cMarker = sPart 
            continue  
         ipos = sPart.find('}')
         if ipos >= 0:           
            sKey = '$' + cMarker + sPart[0:ipos] # try as identifier 
            if dictInlineAnn.has_key(sKey):
               sModified += '{' + dictInlineAnn[sKey] + '}'
            elif self.dictClientAnn.has_key(sKey):
               sModified += '{' + self.dictClientAnn[sKey] + '}'
            else: # sKey failes as identifier, therefore leave unchanged
               sModified += '{' + sKey + '}'
            sModified += sPart[ipos+1:] # add substr occurring left of '}'
         else:
           sMsg = "replace_by_marker: missing '}' after '%s{$*'" % \
                  sModified
           self.lstErrors.append(sMsg)
           return sModified           
      return sModified

   #===================================================================#
   # ELIMINATE any remaining white space                               #
   #===================================================================#
   """------------------------------------------------------------------
      eliminate_backslashes_and_white_space: eliminate backslashes and
                                             and following white space  
                                                                
      return: string with modified notation      
   """
   def eliminate_backslashes_and_white_space(self,sCurrNotation):

      sModified = sCurrNotation
      lstParts  = sCurrNotation.split('\\')
      if len(lstParts) > 1:
         sModified = lstParts[0].strip()
         for sPart in lstParts[1:]:
            sModified += sPart.strip()
               
      # check for blanks
      ipos = sModified.find(' ')
      if ipos > -1:
         sMsg  = 'eliminate_backslashes_and_white_space: '
         sMsg += "unexpected blank at position %d \n in '%s'"\
                 % (ipos,sModified)
         self.lstErrors.append(sMsg)

      return sModified
      
   #===================================================================#
   # REPLACE client aliases                                            #
   #===================================================================#
   """------------------------------------------------------------------
      replace_client_ali: replace client aliases that are aliases for
                          notation components.
      Replacement code should be found in self.dictClientAli;
      otherwise no replacement.
                                                                
      return: string with modified notation      
   """
   def replace_client_ali(self,sCurrNotation):

      # split into subnotations, i.e. dot-separated parts
      lstSub = self.split_into_parts(sCurrNotation,'.')
      if len(self.lstErrors) > 0:
         sMsg  = 'replace_client_ali: split_into_sub unsuccessful' 
         self.lstErrors.append(sMsg)              
         return sCurrNotation

      sModified = ''
      isub = 0
      for sSub in lstSub:
         isub += 1
         if isub > 1:
            sModified += '.'
         ipos = sSub.find('{$')
         if ipos == 0:
            if len(sSub) < 4:
              sMsg  = 'replace_client_ali: %d. subnotation ' % isub
              sMsg += 'too short for client alias'
              self.lstErrors.append(sMsg)                  
              return sCurrNotation

            # find position of first '}' which ends alias but may be followed
            # by component-anchored annotations and/or multiplier
            iposEnd = sSub.find('}')
            if iposEnd < 0:
               sMsg  = 'replace_client_ali: %d. subnotation: ' % isub
               sMsg += "missing '}'"
               self.lstErrors.append(sMsg)                  
               return sCurrNotation

            sAli  = sSub[1:iposEnd]
            sTail = sSub[iposEnd+1:]
      
            if self.dictClientAli.has_key(sAli): # replacement code found
               sModified += self.dictClientAli[sAli] + sTail
            else:
               sModified += sSub  # put back client alias, since no code 
         else:
            sModified += sSub  # keep this no-client-alias subnotation
            
      return sModified
      
   #===================================================================#
   # REPLACE built-in aliases                                          #
   #===================================================================#
   """------------------------------------------------------------------
      replace_builtin_ali: replace built-in aliases that are aliases for
                           notation components.
      Replacement code should be found in modules in /aliases or ;
      otherwise no replacement.
                                                                
      return: string with modified notation      
   """
   def replace_builtin_ali(self,sCurrNotation):

      # split into dot-separated parts
      lstSub = self.split_into_parts(sCurrNotation,'.')
      if len(self.lstErrors) > 0:
         sMsg  = 'replace_builtin_ali: split_into_sub unsuccessful' 
         self.lstErrors.append(sMsg)              
         return sCurrNotation

      sModified = ''
      icompnt = 0
      for sSub in lstSub:
         
         icompnt += 1
         if icompnt > 1:
            sModified += '.'
         
         if len(sSub) < 4:  # ignore one-char "aliases"
            sModified += sSub
            continue

         sCompnt = None
         if sSub[0]=='{' and sSub[1].isalpha():
            iposEnd = sSub.find('}')
            if iposEnd > 2:
               sAlias = sSub[1:iposEnd]
               sCompnt = self.oDataFace.compnt_notation_by_alias(sAlias)
               
         if sCompnt != None:
            sModified += sCompnt + sSub[iposEnd+1:]
         else:
            sModified += sSub

      return sModified


   #===================================================================#
   # MAKE work notation                                          #
   #===================================================================#
   """------------------------------------------------------------------
      make_work_notation: turn current notation into work notation,
                          assuming that the following methods have 
      already been applied to the initial user notation:
      replace_client_ann, eliminate_backslashes_and_white_space,
      replace_client_ali and replace_builtin_ali

      assignments: self.sWorkNotation
                   self.nCompnt  
                   self.lstCompnt  
                   self.lstTypes   
                   self.lstObjSmi (each obj initialized, no parsing yet)
                   self.lstObjSfn (each obj initialized, no parsing yet)

                   self.nCpsParts          
                   self.lstCpsParts   
                   self.lstCpsTypes        
                   self.lstObjSmiCps 
                   self.lstObjSfnCps
                   
      assign: self.sWorkNotation     
      return: self.sWorkNotation
   """
   def make_work_notation(self,sCurrNotation):

      if len(sCurrNotation) < 2:

         # check if upper-case letter for element symbol
         if sCurrNotation.isupper():

            self.sWorkNotation = sCurrNotation
            self.nCompnt = 1  
            self.lstCompnt.append(sCurrNotation)  
            self.lstTypes.append('smi')
            oSmi = csm_annsmi.AnnotatedSmiles(sCurrNotation)
            self.lstObjSmi.append(oSmi)
            self.lstObjSfn.append(None)
            return self.sWorkNotation
         else:
            sMsg  = 'make_work_notation: one-char notation has to' 
            sMsg += ' consist of upper-case letter for organic-set atom'
            self.lstErrors.append(sMsg)                              
            return None

      # Composites and interface-connected structures
      if sCurrNotation[0]=='{' and  sCurrNotation[1]=='/':
         self.nCompnt = 1
         self.sWorkNotation = self.scan_composite_notation(sCurrNotation)
         return self.sWorkNotation

      # apply multiplier and assign remaining variables
      lstCompntTemp = self.split_into_parts(sCurrNotation,'.')
      self.lstCompnt = []
      iCompnt = 0
      self.sWorkNotation = ''          
      for sCompntTemp in lstCompntTemp:
         (sCompnt,sMulti) = self.splitOffRightEndCurly(sCompntTemp)

         nMulti = 1       
         if sMulti != None:
            if sMulti.isdigit():
               nMulti = int(sMulti)
         sSubnotation = sCompntTemp
         if nMulti > 1:
            sSubnotation = sCompnt
         sType = self.evaluate_composite_type(sSubnotation,iCompnt)
         oSmi = oSfn = None         
         for i in range(0,nMulti):
            
            self.lstTypes.append(sType)            

            if cmp(sType,'smi') == 0:
               if i == 0:
                  oSmi = csm_annsmi.AnnotatedSmiles(self.oDataFace,sSubnotation)
                  oSmi.parse()                  
               self.lstObjSmi.append(oSmi)
               self.lstObjSfn.append(None)
               if len(oSmi.msgs_err()) == 0:
                  sMfHill = oSmi.mf_hill_format()
                  self.lstMfHill.append(sMfHill)
               else:
                  self.lstMfHill.append(None)

               self.lstCompntAnn.append(oSmi.caa_entries())
                  
            elif cmp(sType,'sfn') == 0:  
               if i == 0:

                  # split at first '}' to separate sfn from annotations
                  pair = sSubnotation.split('}',1)
                  if len(pair) != 2 or len(pair[0]) < 4:
                     sMsg  = 'make_work_notation: missing "}" in sfn ' 
                     sMsg += '%s (%d. subnotation)' % (sSubnotation,iCompnt)
                     self.lstErrors.append(sMsg)               
                     continue
                  sSfn = pair[0][2:]

                  lstCaa = None
                  # if annotations...
                  if len(pair[1]) > 3:                      
                     lstCaa = self.convert_annseq(pair[1])
                  self.lstCompntAnn.append(lstCaa)
                  
                  oSfn = csm_sfn.StoichFormNotation(self.oDataFace,sSfn)    
                  oSfn.parse()                 
               self.lstObjSmi.append(None)
               self.lstObjSfn.append(oSfn)             
               if len(oSfn.msgs_err()) == 0:
                  sMfHill = oSfn.mf_hill_format()
                  self.lstMfHill.append(sMfHill)
               else:
                  self.lstMfHill.append(None)
            elif cmp(sType,'ali') == 0:
               self.lstObjSmi.append(None)
               self.lstObjSfn.append(None)                           
               self.lstMfHill.append(None)
               self.lstCompntAnn.append(None)
            else:
               self.lstCompntAnn.append(None)


            if iCompnt > 0:
               self.sWorkNotation += '.'            
            self.sWorkNotation += sSubnotation
               
            iCompnt += 1
            
      self.nCompnt = iCompnt
      return self.sWorkNotation


   """------------------------------------------------------------------
      evaluate_component_type: evaluate the component for given
                               component notation
                          
      return: 'smi', 'sfn', 'cps', 'ali', or None
   """
   def evaluate_composite_type(self,sComponent,iCompnt):
      sType = None
      lenCompnt = len(sComponent)
      if lenCompnt < 1:
         sMsg  = 'evaluate_composite_type: %s.' % (iCompnt+1) 
         sMsg += 'component notation has zero length'
         self.lstErrors.append(sMsg)                              
         return sType

      if sComponent[0] != '{':  
         sType = 'smi'
      else:
         if lenCompnt < 3:
            sMsg  = "evaluate_composite_type: %s." % (iCompnt+1)
            sMsg += "component notation, '%s', too short" % sCompnt
            self.lstErrors.append(sMsg)                              
            return None
         if sComponent[1] == '*':
            sType = 'sfn'
         elif sComponent[1] == '/':
            sType = 'cps'
         else:
            sType = 'ali'

      return sType




   """------------------------------------------------------------------
      scan_composite_notation: scan current notation as composite
                               notation and assign variables of
                               composite parts 
                          
      return: sCpsNotation
   """
   def scan_composite_notation(self,sCurrNotation):

      # get annotations
      sEnd = '>'
      sRemainder = sEnd + sCurrNotation
      lstOfPairs = []
      for k in range(len(sCurrNotation)):
         (sRemainder,sStr) = self.splitOffRightEndCurly(sRemainder)
         if cmp(sEnd,sRemainder) == 0:
            sRemainder = '{' + sStr + '}'
            break            
         
         oCurly = csm_curlyann.CurlyAnnotation(self.oDataFace,sStr)
         lstErrorsCurly = oCurly.parse()
         if len(lstErrorsCurly) == 0:
            pair = oCurly.entry()
            lstOfPairs.append(pair)
                            
      self.lstCompntAnn.append(lstOfPairs)
     
      # evaluate composite string
      lstParts = self.split_into_parts(sRemainder[2:-1],'/')
      
      self.nCpsParts = len(lstParts)
      self.lstTypes.append('cps')

      sCpsNotation = '{'
      iPart = 0
      for sPart in lstParts:
         
         iPart += 1
         sCpsNotation += '/'
            
         if sPart[0] == '{':

            if sPart[1] == '*': # SFN 
               self.lstCpsParts.append(sPart)   
               self.lstCpsTypes.append('sfn')
               oSfn = csm_sfn.StoichFormNotation(sCurrNotation)
               self.lstObjSmiCps.append(None)
               self.lstObjSfnCps.append(oSfn)
               sCpsNotation +=  sPart
            else:  # alias
               pass

         else: # (annotated) SMILES notation                   
            self.lstCpsParts.append(sPart)   
            self.lstCpsTypes.append('smi')
            oSmi = csm_annsmi.AnnotatedSmiles(sCurrNotation)
            self.lstObjSmiCps.append(oSmi)
            self.lstObjSfnCps.append(None)
            sCpsNotation += sPart

      sCpsNotation += '}'
       
      return sCpsNotation

   """------------------------------------------------------------------   
      convert a sequence of annotations in Python structure

      EXAMPLE:
         For sAnnSeq = {IMa=Cu}{cr}'
         return:       [['IM', {'a': 'Cu'}], ['cr', {}]]
      
      return: list of pair, where each pair = [AM, {key: val}] 
   """
   def convert_annseq(self, sAnnSeq):
      lstOfPairs = []

      sEnd = '>'
      sRemainder = sEnd + sAnnSeq
      for k in range(len(sAnnSeq)):
         (sRemainder,sAnn) = self.splitOffRightEndCurly(sRemainder)
         oCurly = csm_curlyann.CurlyAnnotation(self.oDataFace,sAnn)
         lstErrorsCurly = oCurly.parse()
         if len(lstErrorsCurly) == 0:
            pair = oCurly.entry()
            lstOfPairs.append(pair)
         
         if cmp(sEnd,sRemainder) == 0:
            break            

      return lstOfPairs

   #===================================================================#
   # HELPER METHODS                                                    #
   #===================================================================#
   """------------------------------------------------------------------   
      split sNotation into parts at cBreak while ignoring cBreak inside
      curly braces
      return: lstParts = list with parts
   """
   def split_into_parts(self,sNotation,cBreak):
      lstParts        = []
      sPart          = ''
      iCurlyDepth = 0 # curly brackets nesting depth
      for cGiv in sNotation:
         if cGiv == cBreak:
            if iCurlyDepth == 0: # split-point reached
               lstParts.append(sPart)
               sPart = ''
               continue
         sPart += cGiv         
         if cGiv == '{':
            iCurlyDepth += 1
         elif cGiv == '}': 
            iCurlyDepth -= 1

      if iCurlyDepth != 0:
         sErr = "split_into_parts: unbalanced curlies in notation"
         self.lstErrors.append(sErr)      
      else:
         lstParts.append(sPart) # append final subnotation

      if len(lstParts) < 1:
         lstParts = [sNotation]

      return lstParts

   """
      splitOffRightEndCurly

      1. Example:  Input:   sStr='ClC=C{Z}Cl{3}'
                   Return:  ('ClC=C{Z}Cl','3')

      2. Example:  Input:   sStr='c1cnccc1'
                   Return:  (None,None)

      3. Example:  Input:   sStr='{NO3(1-)}'
                   Return:  (None,None)                   

      return: (sRemainder,sEndCurly)
   """
   def splitOffRightEndCurly(self,sStr):
      sRemainder  = None # substring of sStr preceding right-end curly
      sEndCurly   = None # content of curly at right end of sStr
      nLen = len(sStr)
      if sStr[nLen-1] == '}':

         # starting at right end of sStr, find nearest '{' to the left,
         # but skip nested pairs of curly braces 
         iPos = nLen-2
         nNestingDepth = 1
         while iPos >= 0:   
            if sStr[iPos] == '}':
               nNestingDepth += 1
            elif sStr[iPos] == '{':
               nNestingDepth -= 1
            if nNestingDepth == 0:
               break
            iPos -= 1
         
         if iPos > 0: # end curly cannot start at position 0
            sRemainder  = sStr[0:iPos]
            sEndCurly   = sStr[iPos+1:nLen-1]
      return (sRemainder,sEndCurly)

   #===================================================================#
   # ACCESS to member variables                                        #
   #===================================================================#

   #===================================================================#
   # MSGS of access violation                                          #
   #===================================================================#
   def icompnt_out_of_range(self,iCompnt,sMeth):
      sMsg  = '%s_compnt(%d): iCompnt out of range \n' % (sMeth,iCompnt)
      sMsg += '   notation has %d components' % self.nCompnt
      self.lstAccFail.append(sMsg)                  
      return None

   def iatom_out_of_range(self,iCompnt,iAtom,sMeth):
      sMsg  = '%s_compnt(%d): iAtom out of range \n' % (sMeth,iAtom)
      sMsg += '   for %d. components' % iCompnt
      self.lstAccFail.append(sMsg)                  
      return None

    
   #===================================================================#
   # ACCESS to SMILES and ANNOTATED SMILES components                  #
   # iuCompnt = iCompnt + 1 (0 < iuCompnt <= self.nCompnt)             #
   #===================================================================#
   def atf0_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'atf0')
   def atoms_aromat_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'aromat')
   def atoms_charge_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'charge')
   def atoms_hcount_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'hcount')         
   def atoms_label_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'label')   
   def atoms_nbors_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'nbors')
   def atoms_nvbors_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'nvbors')
   def atoms_numb_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'numb')
   def atoms_symb_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'symb')   

 
   def aaa_indices_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'aaa_indices')
   def aaa_indices_iu_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'aaa_indices_iu')
   def deloc_charge_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'deloc_charge')   
   def dmat_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'dmat')      
   def msgs_err_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'msgs_err')   
   def numof_atoms_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'numof_atoms')
   def numof_nodes_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'numof_nodes')      
   def numof_rings_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'numof_rings')   
   def rings_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'rings')   
   def rings_iu_compnt(self,iuCompnt):
      return self.csm_compnt(iuCompnt,'rings_iu')   

   def csm_compnt(self,iuCompnt,sAtList):

      if iuCompnt < 1 or iuCompnt > self.nCompnt:
         self.icompnt_out_of_range(iuCompnt,sAtList)
         return None
      else:
         iCompnt = iuCompnt - 1
         if cmp(self.lstTypes[iCompnt],'smi') != 0:
            return None

         #===================================================================#
         # LISTs with data for atom nodes                                    #
         #===================================================================#        
         elif cmp(sAtList,'atf0')==0: 
            return self.lstObjSmi[iCompnt].atf0()      
         elif cmp(sAtList,'aromat')==0:
            return self.lstObjSmi[iCompnt].atoms_aromat()      
         elif cmp(sAtList,'charge')==0:
            return self.lstObjSmi[iCompnt].atoms_charge()     
         elif cmp(sAtList,'hcount')==0:
            return self.lstObjSmi[iCompnt].atoms_hcount()     
         elif cmp(sAtList,'label')==0:
            return self.lstObjSmi[iCompnt].atoms_label()     
         elif cmp(sAtList,'nbors')==0:
            return self.lstObjSmi[iCompnt].atoms_nbors()
         elif cmp(sAtList,'nvbors')==0:
            return self.lstObjSmi[iCompnt].atoms_nvbors()
         elif cmp(sAtList,'numb')==0:
            return self.lstObjSmi[iCompnt].atoms_numb()
         elif cmp(sAtList,'symb')==0:
            return self.lstObjSmi[iCompnt].atoms_symb()

         #===================================================================#
         # other SMILES-component data                                       #
         #===================================================================#   
         elif cmp(sAtList,'aaa_indices')==0:
            return self.lstObjSmi[iCompnt].aaa_indices()
         elif cmp(sAtList,'aaa_indices_iu')==0:
            lst =  self.lstObjSmi[iCompnt].aaa_indices()
            return map(lambda i: i+1, lst)          
         elif cmp(sAtList,'deloc_charge')==0:
            return self.lstObjSmi[iCompnt].deloc_charge()      
         elif cmp(sAtList,'dmat')==0:
            return self.lstObjSmi[iCompnt].dmat()     
         elif cmp(sAtList,'msgs_err')==0:
            return self.lstObjSmi[iCompnt].msgs_err()
         elif cmp(sAtList,'numof_atoms')==0:
            return self.lstObjSmi[iCompnt].numof_atoms() 
         elif cmp(sAtList,'numof_nodes')==0:
            return self.lstObjSmi[iCompnt].numof_nodes()
         elif cmp(sAtList,'numof_rings')==0:
            return self.lstObjSmi[iCompnt].numof_rings()
         elif cmp(sAtList,'rings')==0:
            return self.lstObjSmi[iCompnt].rings()         
         elif cmp(sAtList,'rings_iu')==0:
            return self.lstObjSmi[iCompnt].rings_iu()   

         
         else:
           sMsg  = "atlist_compnt: unknown call  '\n'" % sAtList
           self.lstAccFail.append(sMsg)                  
           return None            

   def aaa_entries_compnt(self,iuCompnt,iuAtom):
      if iuCompnt < 1 or iuCompnt > self.nCompnt:
         self.icompnt_out_of_range(iuCompnt,'aaa_entries_compnt')
         return None
      else:
         nAtoms = self.numof_atoms_compnt(iuCompnt)
         if iuAtom < 1 or iuAtom > nAtoms:
            self.iatom_out_of_range(iuCompnt,iuAtom,'aaa_entries_compnt')
            return None
         else:
            return  (self.lstObjSmi[iuCompnt-1]).aaa_entries(iuAtom)
      

   def atpair_bond_compnt(self,iuCompnt,iuAtom,juAtom):
      if iuCompnt < 1 or iuCompnt > self.nCompnt:
         self.icompnt_out_of_range(iuCompnt,'atpair_bond_compnt')
         return None
      else:
         nAtoms = self.numof_atoms_compnt(iuCompnt)
         if iuAtom < 1 or iuAtom > nAtoms:
            self.iatom_out_of_range(iuCompnt,iuAtom,'atpair_bond_compnt')
            return None
         elif juAtom < 1 or juAtom > nAtoms:
            self.iatom_out_of_range(juCompnt,juAtom,'atpair_bond_compnt')
            return None
         else:
            return  (self.lstObjSmi[iuCompnt-1]).atpair_bond(iuAtom,juAtom) 
   #===================================================================#
   # ACCESS TO COMPONENT ANNOTATIONS                                   #
   #===================================================================#
   def caa_entries_compnt(self,iuCompnt):
      if iuCompnt < 1 or iuCompnt > self.nCompnt:
         self.icompnt_out_of_range(iuCompnt,sAtList)
         return None
      else:
         iCompnt = iuCompnt - 1
         return self.lstCompntAnn[iCompnt]

   #===================================================================#
   # MISCELLANEOUS REQUESTS                                            #
   #===================================================================#
   def msgs_acc(self):              return self.lstAccFail
   def msgs_err(self):              return self.lstErrors
   def numof_components(self):      return self.nCompnt
   def numof_composite_parts(self): return self.nCpsParts
   def type_components(self):       return self.lstTypes
   def type_composite_parts(self):  return self.lstCpsTypes
   def work_notation(self):         return self.sWorkNotation

   def mf_compnt(self,iuCompnt):
      if iuCompnt < 1 or iuCompnt > self.nCompnt:
         self.icompnt_out_of_range(iuCompnt,'mf_compnt')
         return None
      else:
         iCompnt = iuCompnt-1
         sType = self.lstTypes[iCompnt]
         if cmp(sType,'smi') == 0 or cmp(sType,'sfn') == 0:
            return self.lstMfHill[iCompnt]   
         else:
            return None

   def mf_total(self):
      oMfTotal = csm_molform.MolecularFormula(self.oDataFace)
      iCompnt = 0
      for sType in self.lstTypes:
         sMf = None
         if cmp(sType,'smi') == 0:
            oSmi = self.lstObjSmi[iCompnt]
            sMf  = oSmi.mf_linear_notation()
         elif cmp(sType,'sfn') == 0:
            oSfn = self.lstObjSfn[iCompnt]
            sMf  = oSfn.mf_hill_format()
         else:
            return sMf

         # set or add
         if iCompnt == 0:
            oMfTotal.set_molecular_formula(sMf)
            oMfTotal.evaluate()
         else:
            oMfTotal.add_mf(sMf)
         
         iCompnt += 1

      return oMfTotal.hill_format()
