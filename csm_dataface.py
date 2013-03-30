"""
   This file:     csm_dataface.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Python module csm_dataface is intended to be used to
   create a singleton object for look-up of data as nedded
   during parsing and evaluating a CurlySMILES notation.

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
import os
import csm_atoms, csm_aliases

class DataFace:
    
   def __init__(self,sReroute=None):

      self.sCsmpyDir = os.getcwd()
      if sReroute != None:
         self.sCsmpyDir = sReroute
      
      self.oAtoms    = csm_atoms.AtomicData()
      self.oAliases  = csm_aliases.AliasNotations(self.sCsmpyDir)
      self.dictAliases = None # = {sAlias: sPrimAlias,...}
                              #    sAlias = prim. or sec alias
                              
      self.lstStereoDescr = ['D','E','L','R','S','Z']
      self.lstBondSymbols = ['-','=','#','$',':','~']
      
      # group environment annotation markers (GEAMs)
      self.lstGEAMs = ['-R','-X','-Y','=Y','#Y','$Y',':Y','~Y',
                       '-|','~|']
      # molecular detail annotation markers (MDAMs)
      self.lstMDAMs = ['!a','!p','!m','!r','!H','!I']
      
      # operational annotation markers (OPAMs)
      self.lstOPAMs = ['+R','+X','+Y','+L','+n','+r']

      # state and shape annotation markers (SSAMs)
      self.lstSSAMs = ['am', # amorphous
                       'aq', # in aqueous solution
                       'br', # branch or branch-like structure (mesoscale)
                       'cl', # cluster
                       'cp', # colloidal particle
                       'cr', # crystalline
                       'ds', # dissolved (in non-aqueous solvent)
                       'gs', # gaseous
                       'ip', # ion pair
                       'lc', # liquid crystalline
                       'lq', # liquid
                       'ma', # matrix (in composite)
                       'mc', # micelle
                       'mp', # mesoporous
                       'nb', # nanobelt
                       'nc', # nanocrystal
                       'nd', # nanodisk
                       'nf', # nanoflower
                       'nh', # nanowhisker                                   
                       'nk', # nanoflake                       
                       'np', # nanoparticle
                       'nr', # nanorod
                       'ns', # nanostructure
                       'nt', # nanotube
                       'nw', # nanowire
                       'pa', # particle
                       'pc', # poly-crystalline
                       'po', # porous
                       'qd', # quantum dots                       
                       'sc', # single-crystalline
                       'sd', # solid
                       'tf', # thin film
                      ]

      # miscellaneous interest annotation markers (MIAMs)
      self.lstMIAMs = ['AD','IM','SD','@*','@:']

      # annotation dictionary keys
      self.lstAnnDictKeys = ['a','b','c','e','i','n','p',
                             'aa',
                             'all', 'axc','bra','ctr','enz','esa','ful',
                             'hel', 'pep','pha','phn','plm','pro','psy',
                             'sfn', 'spg', 'srf']


      # H-count correction (value double, to have all as integers)
      self.dictHcor = { '-': 2, '=': 4, '#': 6, ':': 3, '$': 8,
                        '-R': 2, '-X': 2, '-Y': 2, '=Y': 4,
                        '#Y': 6, '-|': 2,
                        '+R': 2, '+X': 2, '+Y': 2, '+n': 2, '+r': 2
                      }

      self.dictAliases = None
      self.lstMsgs = []
     
   #===================================================================#
   # ASSESS data availability                                          #
   #===================================================================#
   """
      assess_data_availability: assign self.dictAliases if aliases with
         their corresponding CurlySMILES notations are available
         (if not, place report into self.lstMsgs)
   """
   def assess_data_availability(self):      
      (lstAmbig,dictAliases) = self.oAliases.makeAliasDict()
      if len(lstAmbig) > 0:
         self.lstMsgs += lstAmbig
      else:   
         self.dictAliases = dictAliases
         
   """
       getDataMissReport: get list of lines, each line reporting some
                           data inconsistency or missing data
       return: lstMsgs = [sLine,...] or [], if data availability OK
   """
   def getDataMissReport(self):
      return self.lstMsgs

      
   #===================================================================#
   # QUERY oAtoms                                                      #
   #===================================================================#
   def is_valid_symbol(self,sAtSymb):
      return self.oAtoms.is_valid_symbol(sAtSymb)

   def atnumb_as_str(self,sAtSymb):
      return self.oAtoms.atnumb_as_str(sAtSymb) 
      
   def atnumb_as_int(self,sAtSymb):
      return self.oAtoms.atnumb_as_int(sAtSymb)
 
   def ground_state_electron_configuration(self,sAtSymb):
      return self.oAtoms.getGroundStateElectronConfiguration(sAtSymb)

   def is_valid_one_letter_symbol(self,sAtSymb):
      return self.oAtoms.is_valid_one_letter_symbol(sAtSymb)
      
   def is_in_organic_set(self,sAtSymb):
      return self.oAtoms.is_in_organic_set(sAtSymb)
   
   #===================================================================#
   # QUERY oAliases                                                    #
   #===================================================================#
   """------------------------------------------------------------------
      compnt_notation_by_alias: look up alias in modules 
      
      return: string with component notation that is going to replace
              alias;
              None, if alias not found
   """
   def compnt_notation_by_alias(self, sAlias):
     return  self.oAliases.compnt_notation(sAlias)

   """------------------------------------------------------------------   
       compnt_notation_and_groupid: look up alias and its group-id:
      return:  (ssCompntNotation,sGroupId) or (None,None), if sAlias not found
                sCompntNotation  =  notation for alias replacement
                sGroupId = alias group name such as neutral, cation1p,
                           anion1p,etc.
   """
   def compnt_notation_and_groupid_by_alias(self,sAlias):
       return self.oAliases.compnt_notation_and_groupid(sAlias)

   #===================================================================#
   # QUERY descriptors and annotation markers                          #
   #===================================================================#
   def is_stereodescriptor(self,sStr):
      return sStr in self.lstStereoDescr
   def is_bond_symbol(self,sStr):
      return sStr in self.lstBondSymbols
   def is_GEAM(self,sStr):
      return sStr in self.lstGEAMs
   def is_MDAM(self,sStr):
      return sStr in self.lstMDAMs 
   def is_OPAM(self,sStr):
      return sStr in self.lstOPAMs
   def is_SSAM(self,sStr):
      return sStr in self.lstSSAMs  
   def is_MIAM(self,sStr):
      return sStr in self.lstMIAMs

   def marker_type(self,sStr):
      if self.is_stereodescriptor(sStr):
         return 'SD'
      elif self.is_bond_symbol(sStr):
         return 'BS'
      elif self.is_GEAM(sStr):
         return 'GEAM'
      elif self.is_MDAM(sStr):
         return 'MDAM'
      elif self.is_OPAM(sStr):
         return 'OPAM'
      elif self.is_SSAM(sStr):
         return 'SSAM'
      elif self.is_MIAM(sStr):
         return 'MIAM'      
      else:
         return None

   def is_aaa(self,sStr):
      sMarkerType = self.marker_type(sStr)
      return sMarkerType in ['SD','BS','GEAM','MDAM','OPAM']

   def is_caa(self,sStr):
      sMarkerType = self.marker_type(sStr)
      return sMarkerType in ['SSAM','MIAM']

   def is_open_single_bond(self,sAM):
      return sAM in ['-','-R','-X','-Y','-|',
                     '+R','+X','+Y','+n','+r']

   def is_open_double_bond(self,sAM):
      return sAM in ['=','=Y']

   def is_open_triple_bond(self,sAM):
      return sAM in ['#','#Y']

   def is_open_quadruple_bond(self,sAM):
      return sAM in ['$','$Y']

   def is_open_aromatic_bond(self,sAM):
      return sAM in [':',':Y']

   def is_open_unspecified_bond(self,sAM):
      return sAM in ['~','~|']
      
   #===================================================================#
   # QUERY annotation dictionary keys                                  #
   #===================================================================#
   def is_anndict_key(self,sStr):
      return sStr in self.lstAnnDictKeys
   
   #===================================================================#
   # H-count correction                                                #
   #===================================================================#
   """------------------------------------------------------------------
      hcount_correction:  count of H-atoms that formally 'saturate' open
                          bonds 
      return: zero or formal count for sAM in self.dictHcor
   """
   def hcount_correction(self,sAM):
      if self.dictHcor.has_key(sAM):
         return self.dictHcor[sAM]
      else:
         return 0
