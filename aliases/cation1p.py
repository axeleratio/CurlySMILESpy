"""
   This file:     cation1p.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Class AliasCation1p wraps aliases and their associated SMILES 
   notations for cations with charge 1+.

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

         # ammonium cations
         'NH4(1+)':         '[NH4+]',
         
         'MeNH3(1+)':       'C[NH3+]',
         'Me2NH2(1+)':      'C[NH2+]C',
         'Me3NH(1+)':       'C[NH+](C)C',
         'Me4N(1+)':        'C[N+](C)(C)C',

         'EtNH3(1+)':       'CC[NH3+]',
         'Et2NH2(1+)':      'CC[NH2+]CC',
         'Et3NH(1+)':       'CC[NH+](CC)CC',
         'Et4N(1+)':        'CC[N+](CC)(CC)CC',


         # imidazolium
         'Hmim(1+)':        '[nH]1cn(C)cc1{!re=+}',  # 1-methyl

         # 1,3-disubstituted imidazolium            
         'Amim(1+)':        'C=CCn1cn(C)cc1{!re=+}', # 1-allyl-3-methyl
         'Vmim(1+)':        'C=Cn1cn(C)cc1{!re=+}',  # 1-methyl-3-vinyl       
     
         'mmim(1+)':        'Cn1cn(C)cc1{!re=+1}',

         # 1-methyl-3-alkyl-imidazolium
         'C2mim(1+)':       'CCn1cn(C)cc1{!re=+}',
         'C3mim(1+)':       'CCCn1cn(C)cc1{!re=+}',
         'C4mim(1+)':       'CCCCn1cn(C)cc1{!re=+}',
         'C5mim(1+)':       'CCCCCn1cn(C)cc1{!re=+}',
         'C6mim(1+)':       'CCCCCCn1cn(C)cc1{!re=+}',
         'C7mim(1+)':       'CCCCCCCn1cn(C)cc1{!re=+}',
         'C8mim(1+)':       'CCCCCCCCn1cn(C)cc1{!re=+}',
         'C9mim(1+)':       'CCCCCCCCCn1cn(C)cc1{!re=+}',
         'C10mim(1+)':      'CCCCCCCCCCn1cn(C)cc1{!re=+}',
         'C11mim(1+)':      'CCCCCCCCCCCn1cn(C)cc1{!re=+}',
         'C12mim(1+)':      'CCCCCCCCCCCCn1cn(C)cc1{!re=+}',
         'C13mim(1+)':      'CCCCCCCCCCCCCn1cn(C)cc1{!re=+}',
         'C14mim(1+)':      'CCCCCCCCCCCCCCn1cn(C)cc1{!re=+}',
         'C15mim(1+)':      'CCCCCCCCCCCCCCCn1cn(C)cc1{!re=+}',
         'C16mim(1+)':      'CCCCCCCCCCCCCCCCn1cn(C)cc1{!re=+}',
         'C17mim(1+)':      'CCCCCCCCCCCCCCCCCn1cn(C)cc1{!re=+}',
         'C18mim(1+)':      'CCCCCCCCCCCCCCCCCCn1cn(C)cc1{!re=+}',

         # 1-ethyl-3-alkyl-imidazolium
         'C2eim(1+)':       'CCn1cn(CC)cc1{!re=+}',
         'C3eim(1+)':       'CCCn1cn(CC)cc1{!re=+}',
         'C4eim(1+)':       'CCCCn1cn(CC)cc1{!re=+}',
         'C5eim(1+)':       'CCCCCn1cn(CC)cc1{!re=+}',
         'C6eim(1+)':       'CCCCCCn1cn(CC)cc1{!re=+}',
         'C7eim(1+)':       'CCCCCCCn1cn(CC)cc1{!re=+}',
         'C8eim(1+)':       'CCCCCCCCn1cn(CC)cc1{!re=+}',
         'C9eim(1+)':       'CCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C10eim(1+)':      'CCCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C11eim(1+)':      'CCCCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C12eim(1+)':      'CCCCCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C13eim(1+)':      'CCCCCCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C14eim(1+)':      'CCCCCCCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C15eim(1+)':      'CCCCCCCCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C16eim(1+)':      'CCCCCCCCCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C17eim(1+)':      'CCCCCCCCCCCCCCCCCn1cn(CC)cc1{!re=+}',
         'C18eim(1+)':      'CCCCCCCCCCCCCCCCCCn1cn(CC)cc1{!re=+}',
         
         'C4mmim(1+)':       'CCCCn1c(C)n(C)cc1{!re=+}',

         # nitronium (nitryl)
         'NO2(1+)':         'O=[N]=O{CHe=+1}',

         # phosphonium
         # trihexyl(tetradecyl)phosphonium 
         'C14Hx3P(1+)':     'CCCCCC[P+](CCCCCC)(CCCCCC)CCCCCCCCCCCCCC',

         # N-alkylpyridinium (1-alkylpyridinium)
         '1MePy(1+)':        'Cn1ccccc1{!re=+}',
         '1EtPy(1+)':        'CCn1ccccc1{!re=+}',         
         '1PrPy(1+)':        'CCCn1ccccc1{!re=+}',
         '1BuPy(1+)':        'CCCCn1ccccc1{!re=+}',         
         '1PePy(1+)':        'CCCCCn1ccccc1{!re=+}',
         '1HxPy(1+)':        'CCCCCCn1ccccc1{!re=+}',
         '1HpPy(1+)':        'CCCCCCCn1ccccc1{!re=+}',
         '1OcPy(1+)':        'CCCCCCCCn1ccccc1{!re=+}',

         # N-alkyl-3-methylpyridinium (1-alkyl-3-methylpyridinium)
         '1Me3MePy(1+)':        'Cn1cc(C)ccc1{!re=+}',
         '1Et3MePy(1+)':        'CCn1cc(C)ccc1{!re=+}',         
         '1Pr3MePy(1+)':        'CCCn1cc(C)ccc1{!re=+}',
         '1Bu3MePy(1+)':        'CCCCn1cc(C)ccc1{!re=+}',         
         '1Pe3MePy(1+)':        'CCCCCn1cc(C)ccc1{!re=+}',
         '1Hx3MePy(1+)':        'CCCCCCn1cc(C)ccc1{!re=+}',
         '1Hp3MePy(1+)':        'CCCCCCCn1cc(C)ccc1{!re=+}',
         '1Oc3MePy(1+)':        'CCCCCCCCn1cc(C)ccc1{!re=+}',

         # N-alkyl-4-methylpyridinium (1-alkyl-4-methylpyridinium)
         '1Me4MePy(1+)':        'Cn1ccc(C)cc1{!re=+}',
         '1Et4MePy(1+)':        'CCn1ccc(C)cc1{!re=+}',         
         '1Pr4MePy(1+)':        'CCCn1ccc(C)cc1{!re=+}',
         '1Bu4MePy(1+)':        'CCCCn1ccc(C)cc1{!re=+}',         
         '1Pe4MePy(1+)':        'CCCCCn1ccc(C)cc1{!re=+}',
         '1Hx4MePy(1+)':        'CCCCCCn1ccc(C)cc1{!re=+}',
         '1Hp4MePy(1+)':        'CCCCCCCn1ccc(C)cc1{!re=+}',
         '1Oc4MePy(1+)':        'CCCCCCCCn1ccc(C)cc1{!re=+}',


         # pyrrolidinium
         '1Me1MePyrr(1+)':   'C[N+]1(C)CCCC1',
         '1Et1MePyrr(1+)':   'CC[N+]1(C)CCCC1',
         '1Pr1MePyrr(1+)':   'CCC[N+]1(C)CCCC1',         
         '1Bu1MePyrr(1+)':   'CCCC[N+]1(C)CCCC1',
         '1Pe1MePyrr(1+)':   'CCCCC[N+]1(C)CCCC1',         
         '1Hx1MePyrr(1+)':   'CCCCCC[N+]1(C)CCCC1',
         '1Hp1MePyrr(1+)':   'CCCCCCC[N+]1(C)CCCC1',
         '1Oc1MePyrr(1+)':   'CCCCCCCC[N+]1(C)CCCC1',

         # trihydrogen cation
         'H3(1+)':          '[H]1~[H]~[H]1{!re=+1}',

         # tropylium 
         'C7H7(1+)':        'c1cccccc1{!re=+1}',

         }

      




