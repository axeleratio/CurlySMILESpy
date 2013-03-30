"""
   This file:     code_snippet.py
   Last modified: October 21, 2010
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

if __name__ == '__main__':
   
   oDataFace = csm_dataface.DataFace()
   # if sCsmPath was appended above, then replace previous line with:
   # oDataFace = csm_dataface.DataFace(sCsmPath)
   
   sNotation = "CCCCn1cn(C)cc1{!re=+}.{NTf2(1-)}" 
   oNotation = csm_notation.Notation(oDataFace)
   lstErrors = oNotation.parse(sNotation)
   
   if len(lstErrors) == 0:
      sWorkNotation = oNotation.work_notation()
      sMFTotal  = oNotation.mf_total()
      sMFCation = oNotation.mf_compnt(1)
      print 'sWorkNotation:'
      print '   %s' % sWorkNotation
      print 'sMFTotal:  %s' % sMFTotal

      nComponents = oNotation.numof_components()
      if nComponents > 1:
         sMFAnion  = oNotation.mf_compnt(2)
         lstAnionAtf0 = oNotation.atf0_compnt(2)
         lstAnionDM   = oNotation.dmat_compnt(2)         
         print 'sMFCation: %s' % sMFCation
         print 'sMFAnion:  %s' % sMFAnion
         print ' '
         print 'lstAnionAtf0:'
         print '   %s, %s, %s, %s' % ( lstAnionAtf0[0], lstAnionAtf0[1],\
                                     lstAnionAtf0[2], lstAnionAtf0[3] ) 
         print '   %s, %s, %s, %s,' % ( lstAnionAtf0[4], lstAnionAtf0[5],\
                                     lstAnionAtf0[6], lstAnionAtf0[7] )
         print '   %s, %s, %s, %s,' % ( lstAnionAtf0[8], lstAnionAtf0[9],\
                                     lstAnionAtf0[10], lstAnionAtf0[11] )
         print '   %s, %s, %s' % ( lstAnionAtf0[12], lstAnionAtf0[13],\
                                  lstAnionAtf0[14] )       
         print ' '

         # print upper diagonal of distance matrix row-wise
         print 'lstAnionDM:'
         jcol = oNotation.numof_nodes_compnt(2) 
         lenrow = jcol - 1
         idxAnionDM = 0
         s3Blanks = '   '
         sIndent = s3Blanks
         while lenrow > 0:
            sRow = '%s%d' % (sIndent,lstAnionDM[idxAnionDM])
            idxAnionDM += 1
            for i in range(lenrow-1):            
               sRow += ', %d' %  lstAnionDM[idxAnionDM]
               idxAnionDM += 1
            print sRow
            sIndent += s3Blanks
            lenrow -= 1
   else:
      print  'Found %d errors:' % len(lstErrors)
      for sError in lstErrors:
         print sError

   
 


   
   
