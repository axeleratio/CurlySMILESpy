"""
   This file:     csm_atoms.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Python module csm_atoms is intended to be used to
   look-up atomic data in the context of SMILES and
   CurlySMILES notations.

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
class AtomicData:
    
   def __init__(self):
      self.dictSymbolData      = init_dict_symbol_data()
      self.tplOneLetterSymbols = self.init_tpl_one_letter_symbols()      
      self.tplOrganicSet       = self.init_tpl_organic_set()
      
   #===================================================================#
   # MAKE tuples with special atomic symbols                           #
   #===================================================================#
   """
      init_tpl_one_letter_symbols: initialize and return tuple that
                                   contains all element symbols 
                                   consisting of one letter
   """
   def init_tpl_one_letter_symbols(self):
      tpl = ('H','B','C','N','O','F','P','S','K','V','Y','I','U')
      return tpl

   """
      init_tpl_organic_set: initialize and return tuple with symbols
                            of chemical elements that belong to
                            the organic set (atoms not in the organic
                            set always need to be encoded inside square
                            brackets in SMILES notations)
   """
   def init_tpl_organic_set(self):
      tpl = ('B','C','N','O','F','P','S','Cl','Br','I')
      return tpl

   #===================================================================#
   # ACCESS member objects                                             #
   #===================================================================#  
   def dict_symbol_data(self):
      return self.dictSymbolData

   def tpl_one_Letter_symbols(self):      
      return self.tplOneLetterSymbols

   def tpl_organic_set(self):
      return self.tplOrganicSet

   #===================================================================#
   # QUERY dictSymbolData                                              #
   #===================================================================#
   def is_valid_symbol(self,sAtSymb):
      return self.dictSymbolData.has_key(sAtSymb)

   def atnumb_as_str(self,sAtSymb):
      if sAtSymb == '*':
         return '0'      
      elif self.dictSymbolData.has_key(sAtSymb):
         return self.dictSymbolData[sAtSymb][0]
      else:
         return None
      
   def atnumb_as_int(self,sAtSymb):
      if sAtSymb == '*':
         return 0
      elif self.dictSymbolData.has_key(sAtSymb):
         sAtNumb = self.dictSymbolData[sAtSymb][0]
         return int(sAtNumb)
      else:
         return None

   def ground_state_electron_configuration(self,sAtSymb):
      if self.dictSymbolData.has_key(sAtSymb):
         return self.dictSymbolData[sAtSymb][1]
      else:
         return None 

   #===================================================================#
   # QUERY tplOneLetterSymbols                                         #
   #===================================================================#
   def is_valid_one_letter_symbol(self,sAtSymb):
      return sAtSymb in self.tplOneLetterSymbols

   #===================================================================#
   # QUERY tplOrganicSet                                               #
   #===================================================================#
   def is_in_organic_set(self,sAtSymb):
      return sAtSymb in self.tplOrganicSet


#======================================================================#
# MAKE dictionary with chemical-element data                           #
#======================================================================#
"""
      init_dict_symbol_data: initialize and return dictionary with
                             key   = atomic symbol
                             value = (atomic number,gsec)

         gsec = ground state electron configuration in LaTeX math format
         gsecs are taken from John Emsley: The Elements. Third Edition.
         Oxford University Press, New York, 1998.
"""
def init_dict_symbol_data():
      dictAt = {
         '*':   (   '*', '*'),    # placeholder
         'H':   (   '1', '1s^1' ),
         'He':  (   '2', '1s^2' ),
         'Li':  (   '3', '[He]2s^1' ),        
         'Be':  (   '4', '[He]2s^2' ),          
         'B':   (   '5', '[He]2s^22p^1' ),  
         'C':   (   '6', '[He]2s^22p^2' ),  
         'N':   (   '7', '[He]2s^22p^3' ),  
         'O':   (   '8', '[He]2s^22p^4' ),  
         'F':   (   '9', '[He]2s^22p^5' ),  
         'Ne':  (  '10', '[He]2s^22p^6' ),
         'Na':  (  '11', '[Ne]3s^1' ),          
         'Mg':  (  '12', '[Ne]3s^2' ), 
         'Al':  (  '13', '[Ne]3s^23p^1' ),  
         'Si':  (  '14', '[Ne]3s^23p^2' ),  
         'P':   (  '15', '[Ne]3s^23p^3' ),  
         'S':   (  '16', '[Ne]3s^23p^4' ),  
         'Cl':  (  '17', '[Ne]3s^23p^5' ),  
         'Ar':  (  '18', '[He]3s^23p^6' ),
         'K':   (  '19', '[Ar]4s^1' ),         
         'Ca':  (  '20', '[Ar]4s^2' ),
         'Sc':  (  '21', '[Ar]3d^14s^2' ),
         'Ti':  (  '22', '[Ar]3d^24s^2' ),          
         'V':   (  '23', '[Ar]3d^34s^2' ),   
         'Cr':  (  '24', '[Ar]3d^54s^1' ),
         'Mn':  (  '25', '[Ar]3d^54s^2' ),         
         'Fe':  (  '26', '[Ar]3d^64s^2' ),
         'Co':  (  '27', '[Ar]3d^74s^2' ),
         'Ni':  (  '28', '[Ar]3d^84s^2' ),
         'Cu':  (  '29', '[Ar]3d^{10}4s^1' ),
         'Zn':  (  '30', '[Ar]3d^{10}4s^2' ),
         'Ga':  (  '31', '[Ar]3d^{10}4s^24p^1' ),            
         'Ge':  (  '32', '[Ar]3d^{10}4s^24p^2' ),      
         'As':  (  '33', '[Ar]3d^{10}4s^24p^3' ), 
         'Se':  (  '34', '[Ar]3d^{10}4s^24p^4' ), 
         'Br':  (  '35', '[Ar]3d^{10}4s^24p^5' ),
         'Kr':  (  '36', '[Ar]3d^{10}4s^24p^6' ),
         'Rb':  (  '37', '[Kr]5s^1' ),          
         'Sr':  (  '38', '[Kr]5s^2' ),
         'Y':   (  '39', '[Kr]4d^15s^2' ),         
         'Zr':  (  '40', '[Kr]4d^25s^2' ),          
         'Nb':  (  '41', '[Kr]4d^45s^1' ), 
         'Mo':  (  '42', '[Kr]4d^55s^1' ),
         'Tc':  (  '43', '[Kr]4d^55s^2' ),         
         'Ru':  (  '44', '[Kr]4d^75s^1' ),  
         'Rh':  (  '45', '[Kr]4d^85s^1' ),  
         'Pd':  (  '46', '[Kr]4d^{10}' ),  
         'Ag':  (  '47', '[Kr]4d^{10}5s^1' ),  
         'Cd':  (  '48', '[Kr]4d^{10}5s^2' ),
         'In':  (  '49', '[Kr]4d^{10}5s^25p^1' ),          
         'Sn':  (  '50', '[Kr]4d^{10}5s^25p^2' ),  
         'Sb':  (  '51', '[Kr]4d^{10}5s^25p^3' ),
         'Te':  (  '52', '[Kr]4d^{10}5s^25p^4' ),         
         'I':   (  '53', '[Kr]4d^{10}5s^25p^5' ),          
         'Xe':  (  '54', '[Kr]4d^{10}5s^25p^6' ),    
         'Cs':  (  '55', '[Xe]6s^1' ), 
         'Ba':  (  '56', '[Xe]6s^2' ),
         'La':  (  '57', '[Xe]5d^16s^2' ),         
         'Ce':  (  '58', '[Xe]4f^15d^16s^2' ),
         'Pr':  (  '59', '[Xe]4f^36s^2' ),
         'Nd':  (  '60', '[Xe]4f^46s^2' ),
         'Pm':  (  '61', '[Xe]4f^56s^2' ),
         'Sm':  (  '62', '[Xe]4f^66s^2' ),             
         'Eu':  (  '63', '[Xe]4f^76s^2' ),
         'Gd':  (  '64', '[Xe]4f^75d^16s^2' ),            
         'Tb':  (  '65', '[Xe]4f^96s^2' ),   
         'Dy':  (  '66', '[Xe]4f^{10}6s^2' ),           
         'Ho':  (  '67', '[Xe]4f^{11}6s^2' ),   
         'Er':  (  '68', '[Xe]4f^{12}6s^2' ),
         'Tm':  (  '69', '[Xe]4f^{13}6s^2' ),
         'Yb':  (  '70', '[Xe]4f^{14}6s^2' ),         
         'Lu':  (  '71', '[Xe]4f^{14}5d^16s^2' ),
         'Hf':  (  '72', '[Xe]4f^{14}5d^26s^2' ),
         'Ta':  (  '73', '[Xe]4f^{14}5d^36s^2' ),
         'W':   (  '74', '[Xe]4f^{14}5d^46s^2' ),         
         'Re':  (  '75', '[Xe]4f^{14}5d^56s^2' ),    
         'Os':  (  '76', '[Xe]4f^{14}5d^66s^2' ),   
         'Ir':  (  '77', '[Xe]4f^{14}5d^76s^2' ),    
         'Pt':  (  '78', '[Xe]4f^{14}5d^96s^1' ),
         'Au':  (  '79', '[Xe]4f^{14}5d^{10}6s^1' ),  
         'Hg':  (  '80', '[Xe]4f^{14}5d^{10}6s^2' ),
         'Tl':  (  '81', '[Xe]4f^{14}5d^{10}6s^26p^1' ),         
         'Pb':  (  '82', '[Xe]4f^{14}5d^{10}6s^26p^2' ),
         'Bi':  (  '83', '[Xe]4f^{14}5d^{10}6s^26p^3' ),
         'Po':  (  '84', '[Xe]4f^{14}5d^{10}6s^26p^4' ),         
         'At':  (  '85', '[Xe]4f^{14}5d^{10}6s^26p^5' ),
         'Rn':  (  '86', '[Xe]4f^{14}5d^{10}6s^26p^6' ),
         'Fr':  (  '87', '[Rn]7s^1' ),
         'Ra':  (  '88', '[Rn]7s^2' ),            
         'Ac':  (  '89', '[Rn]6d^17s^2' ),
         'Th':  (  '90', '[Rn]6d^27s^2' ),
         'Pa':  (  '91', '[Rn]5f^26d^17s^2' ),
         'U':   (  '92', '[Rn]5f^36d^17s^2' ),         
         'Np':  (  '93', '[Rn]5f_46d^17s^2' ),         
         'Pu':  (  '94', '[Rn]5f_67s^2' ),  
         'Am':  (  '95', '[Rn]5f^77s^2' ),
         'Cm':  (  '96', '[Rn]5f^76d^17s^2' ),
         'Bk':  (  '97', '[Rn]5f^97s^2' ),
         'Cf':  (  '98', '[Rn]5f^{10}7s^2' ),         
         'Es':  (  '99', '[Rn]5f^{11}7s^2' ),
         'Fm':  ( '100', '[Rn]5f^{12}7s^2' ),
         'Md':  ( '101', '[Rn]5f^{13}7s^2' ),         
         'No':  ( '102', '[Rn]5f^{14}7s^2' ),         
         'Lr':  ( '103', '[Rn]5f^{14}6d^17s^2' ),            
         'Rf':  ( '104', '[Rn]5f^{14}6d^27s^2' ),               
         'Db':  ( '105', '[Rn]5f^{14}6d^37s^2' ),     
         'Sg':  ( '106', '[Rn]5f^{14}6d^47s^2' ),  
         'Bh':  ( '107', '[Rn]5f^{14}6d^57s^2' ),
         'Hs':  ( '108', '[Rn]5f^{14}6d^67s^2' ),
         'Mt':  ( '109', '[Rn]5f^{14}6d^77s^2' )         
      }
      return dictAt
