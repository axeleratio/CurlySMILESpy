"""
   This file:     csm_sfn.py
   Last modified: November 7, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Python module csm_sfn implements a class for managing
   a stoichiometric formula notation.

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
import csm_dataface, csm_molform
import re
class StoichFormNotation:
    
   def __init__(self,oDataFace=None,sSfn=None):

      #
      self.oDataFace = oDataFace
      self.sSfn      = sSfn # re
      self.sSfnChargeStripped = None

      self.lstAtoms = []  # lstPairs = [[sAtSymb,iSubscript],...]
                          # containing atomic-symbol/subscript pairs
                          # of all non-group atoms

      self.lstGroups = [] # lstGroups = [[sGroup,iSubscript],...]
                          # containing group-notation/subscript pairs
                          # of all groups
                         
      self.sCharge = None # charge notation
      self.iCharge = None # charge number
      
      self.oMf = None     # compacted molecular formula: no groups,
                          #   each atomic symbols occurs only once  
      
      self.lstErrors = [] # list of strings with reported error 

   #===================================================================#
   # PARSE SFN                                                         #
   #===================================================================#

   """------------------------------------------------------------------
      parse: parse SFN
             preprocessing: turn user into work notation
             evaluating:

      argument: sUserNotation, to be assigned to self.sUserNotation
                               unless done with self.set_user_notation

      return: error_msgs
   """
   def parse(self):

      # initialize as if no charge notation
      self.sCharge = '0'
      self.iCharge = 0
      self.sSfnChargeStripped = self.sSfn

      # separate charge notation
      # pattern: opening round bracket, followed by one or more digits,
      #    followed by '+' or '-', followed by closing round bracket 
      pat =re.compile('\(([0-9]+)(\+|\-)\)') 
      lstParts = pat.split(self.sSfn)
      
      if len(lstParts) > 2:
         self.sSfnChargeStripped = lstParts[0]
         sNumb = lstParts[1]
         cSign = lstParts[2]
         if not sNumb.isdigit():
            sMsg  = "parse: '%s' invalid charge number" % sNumb
            self.lstErrors.append(sMsg)
            return self.lstErrors
         if not ( cSign == '+' or cSign == '-' ):
            sMsg  = "parse: '%s' invalid charge sign" % cSign
            self.lstErrors.append(sMsg)
            return self.lstErrors                         
         self.sCharge = sNumb + cSign
         self.iCharge = int(sNumb)
         if cSign == '-':
             self.iCharge *= -1


      #
      # pattern: any opening round bracket or a closing round bracket
      #   followed by zero or more digits 
      #sStr = 'Mg3(Pb)Si2(CO8S)1765O5(OH)4'
      pat =re.compile('(\(|\)[0-9]*)') 
      lstParts = pat.split(self.sSfnChargeStripped)
      #print 'csm_sfn.parse: lstParts=',lstParts
      bGroup = 0
      pairGroup = []
      for sPart in lstParts:
         if len(sPart) < 1:
            continue
         if sPart == '(':
            bGroup = 1
            continue
         elif sPart[0] == ')':
            iGroupSubscript = 1
            if len(sPart) > 1:
               iGroupSubscript = int(sPart[1:])
            pairGroup.append(iGroupSubscript)
            self.lstGroups.append(pairGroup)
            pairGroup = []
            bGroup = 0
            continue

         if bGroup == 0:
            self.add_to_lst_atoms(sPart)
         else:
            pairGroup.append(sPart)

      # derive molecular formula by combining atoms and groups 
      self.oMf = csm_molform.MolecularFormula(self.oDataFace)
      if len(self.lstAtoms) > 0:
         self.oMf.set_lst_of_pairs(self.lstAtoms,self.sCharge)
         self.oMf.evaluate()
      if len(self.lstGroups) > 0:
         for pair in self.lstGroups:
            sGroup = pair[0]
            iSubscript = pair[1]
            oMfGroup = csm_molform.MolecularFormula(self.oDataFace,sGroup)
            oMfGroup.evaluate()
            lstOfPairs = oMfGroup.lst_of_pairs()
            lstMulti = []
            if iSubscript > 1:   
               # multiply atom subscript by group subsript for each group atom
               for pair in lstOfPairs:
                  pairNew = pair
                  pairNew[1] = pairNew[1] * iSubscript   
                  lstMulti.append(pairNew)
            else:
               lstMulti = lstOfPairs
               
            self.oMf.add_lst_of_pairs(lstMulti,'0')
         

      return self.lstErrors

   """------------------------------------------------------------------
      add_to_lst_atoms: parse and add atomic symbol/subscript pairs from
                        non-group parts of notation;
                        no validation of atomic symbols here
      updating: self.lstAtoms
   """
   def add_to_lst_atoms(self,sStr):

      oMf = csm_molform.MolecularFormula(self.oDataFace,sStr)
      oMf.evaluate()
      lstOfPairs = oMf.lst_of_pairs()
      lstErrors = oMf.msgs_err()
      if len(lstErrors) > 0:
         for sError in lstErrors:
            sMsg  = "add_to_lst_atoms: %s" % sError
            self.lstErrors.append(sMsg)
         return       
      self.lstAtoms += lstOfPairs
      #print self.lstAtoms
      #print 'groups:',self.lstGroups


   #===================================================================#
   # ACCESS member variables from application                          #
   #===================================================================#
   def linear_notation(self): return self.sSfn
   def linear_notation_charge_stripped(self):
      return self.sSfnChargeStripped      
   def lst_atoms(self):       return self.lstAtoms 
   def lst_groups(self):      return self.lstGroups
   def charge_notation(self): return self.sCharge
   def charge_number(self):   return self.iCharge 
   def msgs_err(self):        return self.lstErrors 

   #===================================================================#
   # DERIVE data for compacted molecular formula                         #
   #===================================================================#
   """------------------------------------------------------------------
      hill_format: sfn converted and compacted into unique formula 
      return: linear notation with molecular formula in hill format 
   """
   def mf_hill_format(self):
      return self.oMf.hill_format()
   
   """------------------------------------------------------------------
      count_atoms:  count the atoms in the formula
      return: nAtoms, the total number of atoms 
   """
   def count_atoms(self):
      return self.oMf.count_atoms()

   """------------------------------------------------------------------
      list_of_pairs: list of atomic-symbol/subscript pairs for compacted
                     formula
      return: lstOfPairs for compacted formula 
   """
   def mf_list_of_pairs(self):
      return self.oMf.list_of_pairs()

   """------------------------------------------------------------------
      list_of_pairs: dict-of-dict for compacted formula
      return: dictOfDict
   """
   def mf_dict_of_dict(self):
      return self.oMf.dict_of_dict()


