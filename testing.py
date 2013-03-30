"""
   This file:     testing.py
   Last modified: November 9, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   'Hello CurlySMILES!'
   This is a hands-on program demonstrating how to include
   CurlySMILES modules into Python application. 

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

# IF code_snippet.py is in a directory that is different from
# the one in which csm_dataface.py and csm_notation.py are in,
# THEN de-comment the following two lines and provide the
# directory path, for example sCsmPath ='/home/me/curlysmiles'
#
# import sys
# sys.path.append(sCsmPath)
import csm_dataface, csm_notation 

def test_notations():

   # list with CurlySMILES notations to be tested
   dictNotations = {
      '[Co+2].[O-]N(=O)=O{2}.O{6}': 
         ['[Co+2].[O-]N(=O)=O.[O-]N(=O)=O.O.O.O.O.O.O',
          9,'smi.smi.smi.smi.smi.smi.smi.smi.smi','CoH12N2O12'],
      '[Co+2].\ [O-]N(=O)=O{2}.\ O{6}':
         ['[Co+2].[O-]N(=O)=O.[O-]N(=O)=O.O.O.O.O.O.O',
          9,'smi.smi.smi.smi.smi.smi.smi.smi.smi','CoH12N2O12'],
      '{*Cr23C6}':
         ['{*Cr23C6}',1,'sfn','C6Cr23'],
      '{*Bi5(4+)}':
         ['{*Bi5(4+)}',1,'sfn','Bi5(4+)'],
      '{*Cu3(CO3)2(OH)2}':
         ['{*Cu3(CO3)2(OH)2}',1,'sfn','C2H2Cu3O8'],
      '{/N{-}=P{+n}(OCCOCCOC)(OCCOCCOC)/{*ZrO2}}':
         ['{/N{-}=P{+n}(OCCOCCOC)(OCCOCCOC)/{*ZrO2}}',
          1,'cps',None], # composites are non-stoichiometric
      '{/[Ag]/{*SiO2}}':
         ['{/[Ag]/{*SiO2}}',1,'cps',None],
      '{/{*Ag}/{*SiO2}}':
         ['{/{*Ag}/{*SiO2}}',1,'cps',None],
      '[*H2]=[*]':
         ['[*H2]=[*]',1,'smi','H2'], # wildcards excluded from mf      
      'O=C(O)C{R}(O)C{R}(O)C(=O)O':
         ['O=C(O)C{R}(O)C{R}(O)C(=O)O',1,'smi','C4H6O6'],
      'O=C(O)C{S}(O)C{S}(O)C(=O)O':
         ['O=C(O)C{S}(O)C{S}(O)C(=O)O',1,'smi','C4H6O6'],
      'O=C(O)C{R}(O)C{S}(O)C(=O)O':
         ['O=C(O)C{R}(O)C{S}(O)C(=O)O',1,'smi','C4H6O6'],  
      'O=CC{D}(O)CO':
         ['O=CC{D}(O)CO',1,'smi','C3H6O3'],
      'O=CC{R}(O)CO':
         ['O=CC{R}(O)CO',1,'smi','C3H6O3'],
      'O=CC{L}(N)CS':
         ['O=CC{L}(N)CS',1,'smi','C3H7NOS'],
      'O=CC{R}(N)CS':
         ['O=CC{R}(N)CS',1,'smi','C3H7NOS'],      
      'O=CC{R}(O)C{R}(O)CO':
         ['O=CC{R}(O)C{R}(O)CO',1,'smi','C4H8O4'],
      'O=C(O)CCCCCCCC=C{E}CCCCCCCC':
         ['O=C(O)CCCCCCCC=C{E}CCCCCCCC',1,'smi','C18H34O2'],
      'O=C(O)CCCCCCCC=C{Z}CCCCCCCC':
         ['O=C(O)CCCCCCCC=C{Z}CCCCCCCC',1,'smi','C18H34O2'],
      'N{-}':
         ['N{-}',1,'smi','H2N'],
      'C{-}#CC#C{-}':
         ['C{-}#CC#C{-}',1,'smi','C4'],
      '[Pt]{~}{~}{~}{~}':
         ['[Pt]{~}{~}{~}{~}',1,'smi','Pt'],
      'C{-Rbra=0}(=O)O{-Rn=4-8}':
         ['C{-Rbra=0}(=O)O{-Rn=4-8}',1,'smi','CO2'],
      'FC{-Yaa=C,Si,Ge}(F)(F)':
         ['FC{-Yaa=C,Si,Ge}(F)(F)',1,'smi','CF3'],      
      'S=C{-Yc=c{:}{:}{-}}':
         ['S=C{-Yc=c{:}{:}{-}}',1,'smi','CHS'],
      '[c]{:}{:}{-Yc=C{-}(=S)}':
         ['[c]{:}{:}{-Yc=C{-}(=S)}',1,'smi','C'],
      '[c]{:}{:}C=S':
         ['[c]{:}{:}C=S',1,'smi','C2HS'],
      'N{=Yc=C{=}(C{-})C{-}}':
         ['N{=Yc=C{=}(C{-})C{-}}',1,'smi','HN'],
      'c1ccccc1{-|c={*CdSe}}':
         ['c1ccccc1{-|c={*CdSe}}',1,'smi','C6H5'],
      'c1ccccc1{-|sfn=CdSe}':
         ['c1ccccc1{-|sfn=CdSe}',1,'smi','C6H5'],      
      'c1n{-R}ccn1{!re=+1}CCC[Si](O)(O{-|sfn=SiO2})O{-|i=-1}':
         ['c1n{-R}ccn1{!re=+1}CCC[Si](O)(O{-|sfn=SiO2})O{-|i=-1}',
          1,'smi','C6H10N2O3Si(1+)'],
      'O=C{!actr=r}(O)C1CC(C{!actr=c}(=O)O)CC(N{!actr=t}(=O)=O)C1':
         ['O=C{!actr=r}(O)C1CC(C{!actr=c}(=O)O)CC(N{!actr=t}(=O)=O)C1',
          1,'smi','C8H11NO6'],
      'Br{!aesa=n}C1CC(C2F{!aesa=s})C(C{!aesa=x}(=O)O)CC21':
         ['Br{!aesa=n}C1CC(C2F{!aesa=s})C(C{!aesa=x}(=O)O)CC21',
          1,'smi','C8H10BrFO2'],
      'O=C(O)C=C{!pi=4;axc=R}C(=O)O':
         ['O=C(O)C=C{!pi=4;axc=R}C(=O)O',
          1,'smi','C4H4O4'],
      'C1C(C(=O)O)CC12CC{!pi=2;axc=S}(C(=O)O)C2':
         ['C1C(C(=O)O)CC12CC{!pi=2;axc=S}(C(=O)O)C2',
          1,'smi','C9H12O4'],      
      'c1cccc1{!re=-}':
         ['c1cccc1{!re=-}',1,'smi','C5H5(1-)'],
      'c1cccccc1{!re=+}':
         ['c1cccccc1{!re=+}',1,'smi','C7H7(1+)'],
      'c1ccccccc1{!re=-2}':
         ['c1ccccccc1{!re=-2}',1,'smi','C8H8(2-)'],
      'CCCCn1cncc1{!re=+}':
         ['CCCCn1cncc1{!re=+}',1,'smi','C7H12N2(1+)'],
      'CCCC[n+]1cn(C)cc1':
         ['CCCC[n+]1cn(C)cc1',1,'smi','C8H15N2(1+)'],
      'CCCCn1c[n+](C)cc1':
         ['CCCCn1c[n+](C)cc1',1,'smi','C8H15N2(1+)'],      
      'O=C(O)c1ccccc1O{!Hi=1}':
         ['O=C(O)c1ccccc1O{!Hi=1}',1,'smi','C7H6O3'],
      'c1ccccc1C(=O)O{!Hi=3#2}.O{!Hi=8#1}C(=O)c1ccccc1':
         ['c1ccccc1C(=O)O{!Hi=3#2}.O{!Hi=8#1}C(=O)c1ccccc1',
          2,'smi.smi','C14H12O4'],
      'c1ccccc1[Si]12OCCN{!Ib=dd;i=7}(CCO2)CCO1':
         ['c1ccccc1[Si]12OCCN{!Ib=dd;i=7}(CCO2)CCO1',1,'smi','C12H17NO3Si'],
      'c1ccccc1[Si]12{!Ib=da;i=11}OCCN(CCO2)CCO1':
         ['c1ccccc1[Si]12{!Ib=da;i=11}OCCN(CCO2)CCO1',1,'smi','C12H17NO3Si'],
      '[Rb+].c1ccccc1[O-]{!Ii=1#1}':
         ['[Rb+].c1ccccc1[O-]{!Ii=1#1}',2,'smi.smi','C6H5ORb'],      
      'N{!I}CCN': ['N{!I}CCN',1,'smi','C2H8N2'],
      'N{!I}CCN{!I}': ['N{!I}CCN{!I}',1,'smi','C2H8N2'],      
      'c1{+Xc=Br{-},Cl{-}}ccc{+R}cc1':
         ['c1{+Xc=Br{-},Cl{-}}ccc{+R}cc1',1,'smi','C6H4'],
      'C1COOCC[N+]1{+Rn=1,2}{+R}':
         ['C1COOCC[N+]1{+Rn=1,2}{+R}',1,'smi','C4H8NO2(1+)'],
      'Cc1ccccc1{+Rc=C{-};p=4,5}':
         ['Cc1ccccc1{+Rc=C{-};p=4,5}',1,'smi','C7H7'],
      'c1ccccc1{+Yc=C{-}n1cncc1;i=1-5}':
         ['c1ccccc1{+Yc=C{-}n1cncc1;i=1-5}',1,'smi','C6H5'],
      '[H]{+Ypep=Thr-Lys-Pro-Arg}':
         ['[H]{+Ypep=Thr-Lys-Pro-Arg}',1,'smi','H'],
      '[H]{+Ypep=TKPR}':
         ['[H]{+Ypep=TKPR}',1,'smi','H'],      
      'c1ccccc1CC(N=[N+]=[N-])C{+Ypep=GKAFVGEI-Met(O)-KS&N{-}}=O':
         ['c1ccccc1CC(N=[N+]=[N-])C{+Ypep=GKAFVGEI-Met(O)-KS&N{-}}=O',1,'smi','C9H8N3O'],
      'c1ccccc1CC(N=[N+]=[N-])C{+Ypep=GKAFVGEI&N{-}C(CCS(=O)C)C(=O)O{-}&KS&N{-}}=O':
         ['c1ccccc1CC(N=[N+]=[N-])C{+Ypep=GKAFVGEI&N{-}C(CCS(=O)C)C(=O)O{-}&KS&N{-}}=O',1,'smi','C9H8N3O'],
      'N{+Yenz=TLL}C(=O)CCc1nnn(c1)CCCN2C(=O)CC(C2=O)CC(C2=O)S{+Ypro=BSA}':
         ['N{+Yenz=TLL}C(=O)CCc1nnn(c1)CCCN2C(=O)CC(C2=O)CC(C2=O)S{+Ypro=BSA}',
          1,'smi','C15H19N5O4S'],      
      '[Co+3].N{!Ii=1#1}{6}':
         ['[Co+3].N{!Ii=1#1}.N{!Ii=1#1}.N{!Ii=1#1}.N{!Ii=1#1}.N{!Ii=1#1}.N{!Ii=1#1}',
          7,'smi.smi.smi.smi.smi.smi.smi','CoH18N6(3+)'],
      '[Co+3]{+Lc=N{6}}':
         ['[Co+3]{+Lc=N{6}}',1,'smi','Co(3+)'],
      '[Hg+2]{+Lc=N#[C-]{!I}{4}}':
         ['[Hg+2]{+Lc=N#[C-]{!I}{4}}',1,'smi','Hg(2+)'],
      '[Pt+2]{+Lc=[Cl-].N{!I}CCN{!I}CCN{!I}}':
         ['[Pt+2]{+Lc=[Cl-].N{!I}CCN{!I}CCN{!I}}',1,'smi','Pt(2+)'],      
      'c1c(C)cc(C)c[ir]1{+Lc=P{!I}(CC)(CC)CC{3}}':
         ['c1c(C)cc(C)c[ir]1{+Lc=P{!I}(CC)(CC)CC{3}}',1,'smi','C7H9Ir'],
      '[W+3]{+Lc=[Cl-]{3}}.[W+3]{+Lc={3}}.[Cl-]{!Ii=1#1,1#2}{3}':
         ['[W+3]{+Lc=[Cl-]{3}}.[W+3]{+Lc={3}}.[Cl-]{!Ii=1#1,1#2}.[Cl-]{!Ii=1#1,1#2}.[Cl-]{!Ii=1#1,1#2}',
          5,'smi.smi.smi.smi.smi','Cl3W2(3+)'],
      '[Ni+]1~[Ni+]~[Ni+]1{!rb=~}{+Lc={C5H5(1-)};i=1,2.[C]{!Ii=1-3}=O{2}}':
      ['[Ni+]1~[Ni+]~[Ni+]1{!rb=~}{+Lc={C5H5(1-)};i=1,2.[C]{!Ii=1-3}=O{2}}',
         1,'smi','Ni3(3+)'],
      'CO{-}CC{+nn=4}OC':
         ['CO{-}CC{+nn=4}OC',1,'smi','C4H8O2'], 
      'N{-}-P{+n}(OCCOCCOC)(OCCOCCOC)':
         ['N{-}-P{+n}(OCCOCCOC)(OCCOCCOC)',1,'smi','C10H24NO6P'],      
      'c1{-}c2OCCOc2c{+n}(s1)':
         ['c1{-}c2OCCOc2c{+n}(s1)',1,'smi','C6H4O2S'],
      'c12OCCOc1c{-}sc2{+ni=1-6}':
         ['c12OCCOc1c{-}sc2{+ni=1-6}',1,'smi','C6H4O2S'],
      'c1c(C(C)(C)C)cc{-}cc1C#C{+rn=6}':
         ['c1c(C(C)(C)C)cc{-}cc1C#C{+rn=6}',1,'smi','C12H12'],
      '[Eu+3]{+Lc=CC(=O)COc1c{-}cc{-R}cc1C{+rn=4}}':
         ['[Eu+3]{+Lc=CC(=O)COc1c{-}cc{-R}cc1C{+rn=4}}',1,'smi','Eu(3+)'],
      'ClCCCCCCC{dsc=ClC(Cl)(Cl)Cl}':
         ['ClCCCCCCC{dsc=ClC(Cl)(Cl)Cl}',1,'smi','C7H15Cl'],
      '[K+]{3}.[O-]P([O-])([O-])=O{aq}':
         ['[K+].[K+].[K+].[O-]P([O-])([O-])=O{aq}',
          4,'smi.smi.smi.smi','K3O4P'],
      'CCCCCCCCCCOc1ccc(cc1)C(=O)Sc2ccc(cc2)CCCCC{lcpha=smC}':
         ['CCCCCCCCCCOc1ccc(cc1)C(=O)Sc2ccc(cc2)CCCCC{lcpha=smC}',
          1,'smi','C28H40O2S'],
      '[C]{crall=diamond}':
        ['[C]{crall=diamond}',1,'smi','C'],
      '[C]{crpsy=cF8}':
        ['[C]{crpsy=cF8}',1,'smi','C'],
      '{*TiO2}':
        ['{*TiO2}',1,'sfn','O2Ti'],
      '{*TiO2}{crphn=anatase}':
        ['{*TiO2}{crphn=anatase}',1,'sfn','O2Ti'],      
      'c1c2c3c4c5c6ccccc6ccc5ccc4ccc3ccc2ccc1{SDhel=P}':
         ['c1c2c3c4c5c6ccccc6ccc5ccc4ccc3ccc2ccc1{SDhel=P}',
          1,'smi','C26H16'],   
      'O=N(=O)c1cccc(C(=O)O)c1-c2c(N(=O)=O)cccc2C(=O)O{SDhel=M}':
         ['O=N(=O)c1cccc(C(=O)O)c1-c2c(N(=O)=O)cccc2C(=O)O{SDhel=M}',
          1,'smi','C14H8N2O8'],
      'O=N(=O)c1cccc(C(=O)O)c1-c2c{!mi=4,12,13;hel=M}(N(=O)=O)cccc2C(=O)O':
         ['O=N(=O)c1cccc(C(=O)O)c1-c2c{!mi=4,12,13;hel=M}(N(=O)=O)cccc2C(=O)O',
          1,'smi','C14H8N2O8'],
      '[2H][Si]([2H])([2H])[2H]{ADsrf=Ni(100)}':
         ['[2H][Si]([2H])([2H])[2H]{ADsrf=Ni(100)}',1,'smi','^2H4Si'],
      '[Na+].[AlH4-]{IMc=[Ti]}':
         ['[Na+].[AlH4-]{IMc=[Ti]}',2,'smi.smi','AlH4Na'],
      '[Ti]{@*c=[Na+].[AlH4-]}':
         ['[Ti]{@*c=[Na+].[AlH4-]}',1,'smi','Ti'],
      '[Ar]{@:ful=C60}':
         ['[Ar]{@:ful=C60}',1,'smi','Ar'],
      'C1CCCCC1{!r$conformation=skew-boat}':
         ['C1CCCCC1{!r$conformation=skew-boat}',1,'smi','C6H12']
   }
   return dictNotations

def evaluate(sNotation,lstExpected): 

   oDataFace = csm_dataface.DataFace()
   # if sCsmPath was appended above, then replace previous line with:
   # oDataFace = csm_dataface.DataFace(sCsmPath)

   sWorkNotationExpected    = lstExpected[0]
   nComponentsExpected      = lstExpected[1]
   sComponentTypesExpected  = lstExpected[2]
   sMFTotalExpected         = lstExpected[3]
   cntDiff = 0
 
   oNotation = csm_notation.Notation(oDataFace)
   lstErrors = oNotation.parse(sNotation)
   if len(lstErrors) == 0:
         sWorkNotationFound = oNotation.work_notation()
         nComponentsFound = oNotation.numof_components()
         lstType = oNotation.type_components()
         sComponentTypesFound = lstType[0]
         for sType in lstType[1:]:
             sComponentTypesFound += '.' + sType   
         sMFTotalFound  = oNotation.mf_total()

         if cmp(sWorkNotationFound, sWorkNotationExpected) != 0:
           print 'Work notations differ for %s' % sNotation
           print '  found:    %s' % sWorkNotationFound
           print '  expected: %s' % sWorkNotationExpected
           cntDiff += 1
         if  nComponentsFound !=  nComponentsExpected:
           print 'Component counts differ for %s' % sNotation
           print '  found:    %s' % nComponentsFound
           print '  expected: %s' % nComponentsExpected               
           cntDiff += 1
         if cmp(sComponentTypesFound, sComponentTypesExpected) != 0:
           print 'Component types differ for %s' % sNotation
           print '  found:    %s' % sComponentTypesFound
           print '  expected: %s' % sComponentTypesExpected
           cntDiff += 1
         if cmp(sMFTotalFound, sMFTotalExpected) != 0:
           print 'MF-totals differ for %s' % sNotation
           print '  found:    %s' % sMFTotalFound
           print '  expected: %s' % sMFTotalExpected
           cntDiff += 1
   else:
      print '%d ERROR(s) for %s' % (len(lstErrors),sNotation)
      for sError in lstErrors:
         print '   %s' % sError
         cntDiff = 9999

   return cntDiff



if __name__ == '__main__':

   # list with CurlySMILES notations to be tested
   dictNotations = test_notations()
   
   cntNotations = 0
   totalDiff = 0
   for sNotation in dictNotations.keys():
      lstExpected = dictNotations[sNotation]
      cntNotations += 1
      cntDiff = evaluate(sNotation,lstExpected)
      if cntDiff > 0:
            totalDiff += 1

   print 'Number of tested notations: %d' % cntNotations 
   print 'Number of notations with found-vs-expected differences: %d'\
         % totalDiff




