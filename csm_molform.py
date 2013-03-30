"""
   This file:     csm_formula.py
   Last modified: October 21, 2010
   Package:       CurlySMILES Version 1.0.1   
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Python module csm_molform implements a class for managing
   a molecular formula. This formula may include isotopically
   labelled atoms and a charge notation at the end.

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
import re
class MolecularFormula:


   def __init__(self,oDataFace=None,sMf=None):
 
      self.oDataFace = oDataFace # reference to object with native data

      """
          A molecular formula, represented with sMf and lstOfPairs allows
          multiple entries of atomic symbols, which may occur in any order.
          A formula based on the Hill convention is obtained with method
          hill_format.
      """
    
      self.sMf        = sMf # molecular formula, linear notation

      self.lstOfPairs = []  # molecular formula as list of pairs,
                            # [sAtSymb,iSubscript]
                             
      self.dictOfDict = {}  # molecular formula as dictionary of
                            # dictionaries:
                            # primary key:   atomic symbol (label-free)
                            # secondary key: labelled atomic symbol
                            # value:         subscript integer ( >0 )

      self.sCharge = None   # charge value notation
      self.iCharge = None   # charge value, whole interger

      self.lstErrors  = []  # list of strings with error messages 

      """
      EXAMPLE-1: (O-D)acetic acid, CH3CO2D
         self.sMf = 'C3H3^2HO2'
         self.lstOfPairs = [['C',3],['H',3],['^2H',1],['O',2]] 
         self.dictOfDict = { 'C': {'C': 3},
                             'H': {'H': 3, '^2H':1},
                             'O': {'O': 2}
                           }
         self.sCharge = '0'                  
         self.iCharge = 0

      EXAMPLE-2: disulfide(2-)
         self.sMf = 'S2(2-)'
         self.lstOfPairs = [ ['S',2] ] 
         self.dictOfDict = { 'S': {'S': 2} }
         self.sCharge = '2-'
         self.iCharge = -2
      """

   #===================================================================#
   # ASSIGNMENT of member variable by application                      #
   #===================================================================#
   """
      set_molecular_formula: assign or reassign self.sMf and
                             initialize all the other member variables
   """
   def set_molecular_formula(self,sMf):
      self.sMf = sMf
      self.lstOfPairs = []   
      self.dictOfDict = {} 
      self.sCharge = '0'  
      self.iCharge = 0   
      self.lstErrors  = []        

   """
      set_list_of_pairs: assign or reassign self.lstOfPairs,
                         self.sCharge and initialize all the
                         other member variables
   """
   def set_lst_of_pairs(self,lstOfPairs,sCharge):
      self.sMf = None
      self.lstOfPairs = lstOfPairs   
      self.dictOfDict = {} 
      self.sCharge = sCharge
      self.iCharge = self.charge_value(sCharge)  
      self.lstErrors  = []        

   """
      set_dict_of_dict: assign or reassign self.dictOfDict,
                        self.sCharge and initialize all the other
                        member variables
   """
   def set_dict_of_dict(self,dictOfDict,sCharge):
      self.sMf = None
      self.lstOfPairs = []
      self.dictOfDict = dictOfDict 
      self.sCharge = sCharge
      self.iCharge = self.charge_value(sCharge)   
      self.lstErrors  = []      

   """
      charge_value: derive charge number (a whole number) from
                    charge notation of type '0', 'n+' or 'n-',
                    where n is a positive interger
      return: iNumber with whole number charge value               
   """
   def charge_value(self,sChargeNotation):
      if sChargeNotation == '0':
         return 0
      else:
         if len(sChargeNotation) < 2:
            sMsg  = "charge_value: '%s' too short, " % sChargeNotation
            sMsg += "minimum length for charge notation is 2 "
            self.lstErrors.append(sMsg)            
            return None
         cSign = sChargeNotation[-1]
         if not ( cSign == '+' or cSign == '-' ):
            sMsg  = "charge_value: '%s' is invalid charge sign" % cSign
            self.lstErrors.append(sMsg)            
            return None
         sNumb = sChargeNotation[0:-1]
         if not sNumb.isdigit():
            sMsg  = "charge_value: '%s' is not all-digits" % sNumb
            self.lstErrors.append(sMsg)            
            return None
         iNumber = int(sNumb)
         if cSign == '-':
            iNumber *= -1
         return iNumber

   #===================================================================#
   # EVALUATE                                                          #
   #===================================================================#
   """------------------------------------------------------------------
      evaluate: evaluate molecular formula depending on entered format,
                assign other member variables
      return: self.lstErrors          
   """
   def evaluate(self):
      if self.sMf != None:  
         self.entry_point_mf()
      elif len(self.lstOfPairs) > 0:
         self.entry_point_lst_of_pairs()
      elif len(self.dictOfDict) > 0:
         self.entry_point_dict_of_dict()
      else:  # nothing to begin with
         pass

      return self.lstErrors

   """------------------------------------------------------------------
      evaluate: evaluate molecular formula from linear notation
   """
   def entry_point_mf(self):

      # separate charge from atomic symbol sequence
      pair = self.sMf.split('(')
      sMfPure = self.sMf
      if len(pair) == 2:
         sMfPure = pair[0]
         if len(pair[1]) < 3:
            sMsg  = 'entry_point_mf: charge notation too short'
            self.lstErrors.append(sMsg)
            return                      
         if pair[1][-1] != ')':
            sMsg  = 'entry_point_mf: missing closing round bracket '
            sMsg += ' for charge notation'
            self.lstErrors.append(sMsg)
            return                                   
         self.sCharge = pair[1][0:-1]
         cSign = self.sCharge[-1]
         sNumb = self.sCharge[0:-1]
         if not sNumb.isdigit():
            sMsg  = "entry_point_mf: '%s' invalid charge number" % sNumb
            self.lstErrors.append(sMsg)
            return
         if cSign == '-':
             self.iCharge = (-1) * int(sNumb)
         elif cSign == '+':   
            self.iCharge = int(sNumb)
         else:
            sMsg  = "entry_point_mf: '%s' invalid charge sign" % cSign
            self.lstErrors.append(sMsg)
            return                                                      
      elif  len(pair) == 1:   
         self.sCharge = '0'                  
         self.iCharge = 0
      else:
         sMsg = \
            'entry_point_mf: found more than one opening round bracket'
         self.lstErrors.append(sMsg)
         return

      self.lstOfPairs = self.parse_atsymb_seq(sMfPure)
      self.dictOfDict = self.lst2dict(self.lstOfPairs)   
      
   """------------------------------------------------------------------
       entry_point_lst_of_pairs: evaluate molecular formula from
                                 lst-of-pairs format
   """
   def entry_point_lst_of_pairs(self):
      self.dictOfDict = self.lst2dict(self.lstOfPairs)
      self.sMf = self.hill_format()
      
          
   """------------------------------------------------------------------
      entry_point_dict_of_dict: evaluate molecular formula from
                                dict-of-dict format
   """
   def entry_point_dict_of_dict(self):
       self.sMf = self.hill_format()
       self.lstOfPairs = self.parse_atsymb_seq(self.sMf) 

   """------------------------------------------------------------------
      parse_atsymb_seq: parse atomic symbol sequence 
      return: lstOfPairs
   """
   def parse_atsymb_seq(self,sSeq):
      lstPairs = []

      # break sStr at each upper-case letter, but retain them
      pat =re.compile('([A-Z]|\^)') 
      lstParts = pat.split(sSeq)
      nParts = len(lstParts)
      
      sAtSymb = ''
      iSubscript = 1
      bExpectLabel = 0
      iParts = 0
      for sPart in lstParts:
         iParts += 1
         if len(sPart) < 1:
            continue

         if sPart[0] == '^':         # isotope label
            sAtSymb += sPart
            bExpectLabel = 1
         elif bExpectLabel == 1:
            sAtSymb += sPart
            bExpectLabel = 0
         elif sPart[0].isupper():    # atomic symbol, first letter
            sAtSymb += sPart
            # look ahead
            if iParts < nParts:
               sAhead = lstParts[iParts]
               if len(sAhead) < 1:
                  iSubscript = 1      # reconfirm 1 for subscript
               elif sAhead.isdigit(): # subscript after one-letter symbol
                  iSubscript = int(sAhead)
               elif len(sAhead) > 0:
                  if sAhead[0].islower(): # atomic symbol, second letter
                     sAtSymb += sAhead[0]
                     if len(sAhead[1:]) > 0: # subscript after two-letter
                                             # symbol
                        if sAhead[1:].isdigit():
                           iSubscript = int(sAhead[1:])
                     else:                   # if no explicit subscript
                        iSubscript = 1       # after two-letter symbol
                     
            pair = [sAtSymb,iSubscript]
            lstPairs.append(pair)
            sAtSymb = ''
            iSubsript = 1
            bExpectLabel = 0
               
      return lstPairs
      

   """------------------------------------------------------------------
      lst2Dict: derive dictOfDict from lstOfPairs 
      return: dictOfDict
   """
   def lst2dict(self,lstOfPairs):
      dictOfDict = {}
      
      for pair in lstOfPairs:
         sAtSymb = pair[0]
         [sLabel,sAtSymbPure] = self.split_labelled_atsymb(sAtSymb)
         iSubscript = pair[1]

         if dictOfDict.has_key(sAtSymbPure):
            subdict = dictOfDict[sAtSymbPure]
            if subdict.has_key(sAtSymb):
               subdict[sAtSymb] =  subdict[sAtSymb] + iSubscript
            else:
               subdict[sAtSymb] = iSubscript
         else:
            dictOfDict[sAtSymbPure] = {sAtSymb: iSubscript}


      return dictOfDict



   def split_labelled_atsymb(self,sAtSymb):
      sLabel = ''
      sSymb  = ''
      bDoSymb = 0
      for cChar in sAtSymb:
         if cChar.isupper():
            bDoSymb = 1
            sSymb += cChar
         elif bDoSymb == 1:   
            sSymb += cChar
         else:
            sLabel += cChar
      return [sLabel,sSymb]





   #===================================================================#
   # ADD                                                               #
   #===================================================================#
   """------------------------------------------------------------------
      add_mf: formally add sMf to current formula by also updating
              lstOfPairs, dictOfDict and charge variables
      return: new mf        
   """
   def add_mf(self,sMf):
      
      # separate charge from atomic symbol sequence
      pair = sMf.split('(')
      sMfPure = sMf
      if len(pair) == 2:
         sMfPure = pair[0]
         if len(pair[1]) < 3:
            sMsg  = 'add_mf: charge notation too short'
            self.lstErrors.append(sMsg)
            return                      
         if pair[1][-1] != ')':
            sMsg  = 'add_mf: missing closing round bracket '
            sMsg += ' for charge notation'
            self.lstErrors.append(sMsg)
            return                                   
         sCharge = pair[1][0:-1]
         cSign = sCharge[-1]
         sNumb = sCharge[0:-1]
         if not sNumb.isdigit():
            sMsg  = "add_mf: '%s' invalid charge number" % sNumb
            self.lstErrors.append(sMsg)
            return
         if cSign == '-':
            iCharge = (-1) * int(sNumb)
         elif cSign == '+':   
            iCharge = int(sNumb)
         else:
            sMsg  = "add_mf: '%s' invalid charge sign" % cSign
            self.lstErrors.append(sMsg)
            return                                                      
      elif  len(pair) == 1:   
         sCharge = '0'                  
         iCharge = 0
      else:
         sMsg = \
            'add_mf: found more than one opening round bracket'
         self.lstErrors.append(sMsg)
         return None

      # update charge
      self.iCharge += iCharge
      if self.iCharge != 0:
         self.sCharge = str(abs(self.iCharge))
         if self.iCharge > 0:
            self.sCharge += '+'
         else:
            self.sCharge += '-'        
      else:
         self.sCharge = '0'

      # make lst of pairs and dict of dict for formula to be added
      lstOfPairs = self.parse_atsymb_seq(sMfPure)
      dictOfDict = self.lst2dict(lstOfPairs)


      # add dictOfDict together
      self.dictOfDict = self.dict_plus_dict(self.dictOfDict, dictOfDict)

      # update remaining members
      self.entry_point_dict_of_dict()

      return self.sMf


   """------------------------------------------------------------------
      add_lst_of_pairs: formally add lst_of_pairs to current formula by
                        updating all member variables accordingly
      return: new sMf        
   """
   def add_lst_of_pairs(self,lstOfPairs,sCharge):
   
      # uppdate charge
      self.iCharge += self.charge_value(sCharge)
      if self.iCharge != 0:
         self.sCharge = str(abs(self.iCharge))
         if self.iCharge > 0:
            self.sCharge += '+'
         else:
            self.sCharge += '-'        
      else:
         self.sCharge = '0'

      # make dict of dict for list of pairs to be added
      dictOfDict = self.lst2dict(lstOfPairs)

      # add dictOfDict together
      self.dictOfDict = self.dict_plus_dict(self.dictOfDict, dictOfDict)

      # formula, representing sum 
      self.sMf = self.hill_format()

      # list of pairs, representing sum
      pair = self.sMf.split('(') # cut-off charge notation
      self.lstOfPairs = self.parse_atsymb_seq(pair[0])
      
      return self.sMf


   """------------------------------------------------------------------
      add_dict_of_dict: formally add dict_of_dict to current formula by
                        updating all member variables accordingly
      return: new sMf        
   """
   def add_dict_of_dict(self,dictOfDict,sCharge):
   
      # uppdate charge
      self.iCharge += self.charge_value(sCharge)
      if self.iCharge != 0:
         self.sCharge = str(abs(self.iCharge))
         if self.iCharge > 0:
            self.sCharge += '+'
         else:
            self.sCharge += '-'        
      else:
         self.sCharge = '0'

      # add dictOfDict together
      self.dictOfDict = self.dict_plus_dict(self.dictOfDict, dictOfDict)

      # formula, representing sum 
      self.sMf = self.hill_format()

      # list of pairs, representing sum
      pair = self.sMf.split('(') # cut-off charge notation
      self.lstOfPairs = self.parse_atsymb_seq(pair[0])
      
      return self.sMf

   
   """------------------------------------------------------------------
      dict_plus_dict: add (merge) two dictOfDict structures
      return: new dictOfDict from combination of dictA and dictB     
   """
   def dict_plus_dict(self,dictA,dictB):
      dictSum = {}     

      # include/add atomic symbols from A and those that occur in A and B
      for sAtSymb in dictA.keys():
         subdictA = dictA[sAtSymb]
         if dictB.has_key(sAtSymb):
            subdictB = dictB[sAtSymb]
            dictSum[sAtSymb] = \
               self.subdict_plus_subdict(subdictA,subdictB) 
         else:
            dictSum[sAtSymb] = subdictA

      # include atomic symbols that occur only in B
      for sAtSymb in dictB.keys():
         if not dictSum.has_key(sAtSymb):         
            subdictB = dictB[sAtSymb]
            dictSum[sAtSymb] = subdictB

      return dictSum

   """------------------------------------------------------------------
      subdict_plus_subdict: add two subdictionaries by adding subscripts
                            in correspondance to atomic symbols

      EXAMPLE:  subdictA   = {'^2H': 5}
                subdictB   = {'H': 3, '^2H': 2}
                subdictSum = {'H': 3, '^2H': 7}
      
      return: new subdictionary from combination of subdictA and subdictB     
   """
   def subdict_plus_subdict(self,subdictA,subdictB):

      subdictSum = {}

      # include/add subscripts from A and those that occur in A and B
      for sSpecific in subdictA.keys():
         iSubscrA = subdictA[sSpecific]
         if subdictB.has_key(sSpecific):
            iSubscrB = subdictB[sSpecific]
            subdictSum[sSpecific] = iSubscrA + iSubscrB
         else:
            subdictSum[sSpecific] = iSubscrA

      # include subscripts that occur only in B
      for sSpecific in subdictB.keys():
         if not subdictSum.has_key(sSpecific):
            iSubscrB = subdictB[sSpecific]
            subdictSum[sSpecific] = iSubscrB

      return subdictSum

   #===================================================================#
   # VALIDATION of atomic symbols                                      #
   #===================================================================#

   """------------------------------------------------------------------
      verify_atsymb: verify atomic symbols in self.dictOfDict     
      return: lstUnknown, list with unknown symbols   
   """
   def verify_atsymb(self):
     lstUnknown = []
     for sAtSymb in self.dictOfDict.keys():
        if not self.oDataFace.is_valid_symbol(sAtSymb):
           lstUnknown.append(sAtSymb)
     return lstUnknown


   #===================================================================#
   # ACCESS member variables from application                          #
   #===================================================================#
   def linear_notation(self): return self.sMf   
   def lst_of_pairs(self):    return self.lstOfPairs 
   def dict_of_dict(self):    return self.dictOfDict
   def charge_notation(self): return self.sCharge
   def charge_number(self):   return self.iCharge 
   def msgs_err(self):        return self.lstErrors 

   #===================================================================#
   # COUNT atoms                                                       #
   #===================================================================#
   """------------------------------------------------------------------
      count_atoms:  count the atoms in the formula
      return: nAtoms, the total number of atoms 
   """
   def count_atoms(self):
      nAtoms = 0
      for pair in self.lstOfPairs:
         nAtoms += pair[1]
      return nAtoms   

   #===================================================================#
   # GENERATE molecular formula in Hill format                         #
   #===================================================================#

   """------------------------------------------------------------------
      hill_format:  generate a molecular formula according to Hill
                    convention, whereafter a formula is written with C
      first, H second, and than all other elements in alphabetical order
      of their atomic symbols. Symbols of the same element are in the
      order: unlabelled first, followed by isotopically labelled symbols,
      starting with lowest label and continuing with increasing label
      number. 
      return: sMfHill
   """
   def hill_format(self):
      sMfHill = ''
      
      lstAtSymb = self.dictOfDict.keys()
      lstAtSymb.sort()

      # symbols for carbon first
      iposC = -1
      ipos = 0
      for sAtSymb in lstAtSymb:
         if sAtSymb == 'C': iposC = ipos 
         ipos += 1
      if iposC > -1:
         subdict = self.dictOfDict['C']
         if len(subdict) > 1:
            sMfHill += self.label_format(subdict) 
         else:
            lst = subdict.keys()
            sSpecific = lst[0]   # the only one, may differ from sAtSymb  
            sMfHill += sSpecific # 'C' or labelled C
            iSubscript = subdict[sSpecific]
            if iSubscript > 1:
               sMfHill += str(iSubscript)      
         del lstAtSymb[iposC]

         # symbols for hydrogen following carbon
         iposH = -1
         ipos = 0
         for sAtSymb in lstAtSymb:
            if sAtSymb == 'H': iposH = ipos 
            ipos += 1
         if iposH > -1:
            subdict = self.dictOfDict['H']
            if len(subdict) > 1:
               sMfHill += self.label_format(subdict) 
            else:
               lst = subdict.keys()
               sSpecific = lst[0]   # the only one, may differ from sAtSymb  
               sMfHill += sSpecific # 'H' or labelled H
               iSubscript = subdict[sSpecific]
               if iSubscript > 1:
                  sMfHill += str(iSubscript)                        
            del lstAtSymb[iposH]

      # all symbols for elements that are not carbon or hydrogen 
      for sAtSymb in lstAtSymb:
         subdict = self.dictOfDict[sAtSymb]
         
         if len(subdict) > 1:
            sMfHill += self.label_format(subdict) 
         else:
            lst = subdict.keys()
            sSpecific = lst[0]   # the only one, may differ from sAtSymb
            sMfHill += sSpecific # bare atomic symbol or labelled
            iSubscript = subdict[sSpecific]
            if iSubscript > 1:
               sMfHill += str(iSubscript)

      # append charge notation, if different from zero
      if self.iCharge != 0:
         if self.iCharge > 0:
            sMfHill += '(%d+)' % self.iCharge
         else:
            sMfHill += '(%d-)' % ( (-1) * self.iCharge )

      return sMfHill
   
   """------------------------------------------------------------------
      label_format: order symbols of the same element: unlabelled first,
                    followed by isotopically labelled symbols and
                    generate linear notation accordingly.

      EXAMPLE: dict = { '^18O' : 2, 'O': 5, '^16O': 1 }
               sSubseq = 'O5^16O^18O2'
                    
      return: sSubseq, linear notation from symbols in dict
   """
   def label_format(self,dict):
      sSubseq = ''

      # make dictionary with key=integer-lables, value=sAtSymb
      dictLabels = {}
      for sAtSymb in dict.keys():
         if sAtSymb[0] == '^':
            [sLabel,sSymbPure] = self.split_labelled_atsymb(sAtSymb)
            dictLabels[int(sLabel[1:])] = sAtSymb
         else:
            dictLabels[0] = sAtSymb # key=0 for unlabelled symbol 

      # sort by label
      lstLabels = dictLabels.keys()
      lstLabels.sort()

      # build linear notation with label-sorted symbols
      for iLabel in lstLabels:
         sAtSymb = dictLabels[iLabel]
         iSubscript = dict[sAtSymb]
         sSubseq += sAtSymb
         if iSubscript > 1:
            sSubseq += str(iSubscript)
         
      return sSubseq
