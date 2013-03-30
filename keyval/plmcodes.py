"""
   This file:     plmcodes.py
   Last modified: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      www.axeleratio.com/axel/axel.htm

   Class PolymerCodes wraps short codes for polymers/plastics,
   excluding biopolymers, with their associated chemical name.
   The short includes acronyms and short brand names.
   To be used as entry in annotation dictionary: {AMplm=...}

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
class PolymerCodes:

   def __init__(self,assignDict=1):      
      self.dictPlm = None
      if assignDict == 1:
         self.initDict()      

   def getSciTerm(self,sAlias):
      if self.dictPlm.has_key(sCode):
         return self.dictPlm[sCode]
      else:
         return None

   def getDict(self):
      return self.dictPlm

   def initDict(self):

      self.dictPlm = {

         'ABR':      'acrylate-butadiene rubber',

         'BR':       'butadiene rubber',

         'CA':       'cellulose acetate',
         'CAB':      'cellulose acetate butyrate',         
         'CAP':      'cellulose acetate propionate',

         'HDPE':     'high-density poly(ethylene)',
         
         'LLDPE':    'linear low-density poly(ethylene)',
         'LMDPE':    'linear medium-density poly(ethylene)',         
         'LDPE':     'low-density poly(ethylene)',

         'NBR':      'acrylonitrile-butadiene rubber (nitrile rubber)',

         'PA':       'aliphatic polyamide (nylon)',
         'PA_4.2':   'poly(tetramethylene oxalamide)',
         'PA_4.6':   'poly(tetramethylene adipamide)',
         'PA_6':     'poly(epsilon-caprolactam)',
         'PA_6.6':   'poly(hexamethylene adipamide)',
         'PA_6.9':   'poly(hexamethylene azellamide)',         
         'PA_6.10':  'poly(hexamethylene sebacamide)',
         'PA_6.12':  'poly(hexamethylene dodecanoamide)',
         'PA_11':    'poly(11-aminoundecanoic acid)',
         'PA_12':    'poly(omega-laurolactam)',
         'PAA':      'poly(acrylic acid)',
         'PADC':     'poly(allyl diglycol carbonate)',
         'PAN':      'poly(acrylonitrile)',
         'PB':       'poly(1-butene)',
         'PBI':      'poly(benzimidazole)',
         'PBR':      'vinyl pyridine-butadiene copolymer',
         'PC':       'polycarbonate',
         'PCTFE':    'poly(chlorotrifluoroethylene)',
         'PDAP':     'poly(diallyl phthalate)',
         'PE':       'poly(ethylene)',
         'PE-C':     'chlorinated poly(ethylene)',
         'PEO':      'poly(ethylene oxide)',
         'PETP':     'poly(ethylene terephthalate)',
         'PI':       'polyimide',
         'PIB':      'poly(isobutylene)',
         'PMMA':     'poly(methyl methacrylate)',
         'PMP':      'poly(4-methyl-1-pentene)',
         'POM':      'poly(oxymethylene)',
         'PPOX':     'poly(propylene oxide)',
         'PSU':      'polysulfone',
         'PTFE':     'poly(tetrafluoroethylene)',
         'PUR':      'polyurethane',
         'PVAC':     'poly(vinyl acetate)',         
         'PVAL':     'poly(vinyl alcohol)',
         'PVB':      'poly(vinyl butyral)',
         'PVC':      'poly(vinyl chloride)',
         'PVCA':     'vinyl chloride-vinyl acetate copolymer',
         'PVDC':     'poly(vinylidene chloride)',
         'PVDF':     'poly(vinylidene fluoride)',
         'PVF':      'poly(vinyl fluoride)',
         'PVFM':     'poly(viny formal)',
         'PVK':      'poly(N-vinyl carbazole)',
         'PVP':      'poly(N-vinyl pyrrolidone)',
         'PVPP':     'polyvinylpolypyrrolidone', # highly cross-linked PVP        

         'SCR':      'styrene-chloroprene copolymer (rubber)',

         'UF':       'urea-formaldehyde resin',
         'UP':       'unsaturated polyester',
      }

      




