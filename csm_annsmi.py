"""
   This file:     csm_annsmi.py
   Last modofied: October 30, 2010
   Package:       CurlySMILES Version 1.0.1
   Author:        Axel Drefahl
   E-mail:        axeleratio@yahoo.com
   Internet:      http://www.axeleratio.com/csm/proj/main.htm

   Python module csm_annsmi implements a class for managing
   a CurlySMILES subnotation (component notation) of type 'smi'.
   The subnotation string is expected to contain a SMILES (with
   no dots) or a CurlySMILES component encoding.

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
import csm_dataface, csm_curlyann, csm_molform

class AnnotatedSmiles:
    
   def __init__(self,oDataFace=None,sAnnSmi=None):

      self.oDataFace = oDataFace

      # CurlySMILES notation (component of a work notation)
      self.sAnnSmi   = sAnnSmi 

      # variables assigned during preparsing
      self.lstTokens = None
      self.lstTokTyp = None
      self.lstDepth  = None
      self.maxDepth  = None

      # variables assigned during parsing
      #   (1) counts
      self.nNodes = None # number of nodes in molecular graph, includes all atoms
                         # in notation, except those H-atom that follow another
                         # atomic symbol inside square-bracket notation
      self.nAtoms = None # nNodes + nHterm
      self.nHterm = None # number non-node H-atoms adjacent to node atom
      self.nAhold = None # number of atomic symbol placeholders      
      self.nRings = None # number of rings
      self.nLocal = None # localized charge count over all atoms      
      self.nDeloc = None # delocalized charge count

      #   (2) attributes of nodes in molecular graph
      #       (node = atom with atomic symbol occurring in notation; notice that
      #               hydrogen atoms can be suppressed)
      self.lstAtSymb    = None # list of str with atomic symbols
      self.lstAtNumb    = None # list of int with atomic number
      self.lstAtLabel   = None # list of int with isotope label
      self.lstAtCharge  = None # list of int with charge value
      self.lstAtDepth   = None # list of int with nesting depth    
      self.lstAromat    = None # list of boolean (1=aromatic-ring atom, 0 otherwise)
      self.lstAaaEntr   = None # list of list with CurlyAnnotation objects
                               # for atom-anchored annotations (AAA)
                               # (empty list for non-annotated atoms)
      self.lstHcount    = None # list of int with number of adjacent non-node H-atoms
      self.lstNbors     = None # list of int with number of neighbor nodes      
      self.lstNvbors    = None # list of int with number of valence bonds to
                               # neighbor nodes    
      self.lstRingMap   = None # list of indices indicating rings of an atom
                               # is a member 
      self.lstCaaEntr   = None # list of list with CurlyAnnotation objects
                               # for component-anchored annotations (CAA)
      self.lstAtf0      = None # atom fragments (level 0: central atom plus
                               # adjacent non-node H-atoms plus dangling bonds)

      #   (3) connectivity
      self.dictNbors   = None # dictionary with neigbors:
                              #    {idxAt:[i0,i1,...],...}
                              #    idxAt   = index of atom;
                              #    i0,i1,..= indexes of its neighbors
      self.dictBonds   = None # dictionary with bonds to neigbors:
                              #    {idxAt:[s0,s1,...],...}
                              #    idxAt   = index of atom;
                              #    s0,a1,..= bond symbols                 
      self.dictPairs   = None # dictionary with bond between bound atom pair
                              #    {'i,j': sBond,...}
                              #     i,j  = atom indexes for pair; i<j
                              #    sBond = bond symbol

      #   (4) rings
      self.lstRingLength = None
      self.lstRingAromat = None
      self.lstRingCharge = None
      self.lstRings      = None

      #   (5) topological matrices
      self.lstIdxMat  = None
      self.lstAdjMat  = None
      self.lstDistMat = None

      #   (6) molecular formula
      self.oMf  = None

      #   (7) status information 
      self.lstErrors = []
      
   #===================================================================#
   # INITIALIZE                                                        #
   #===================================================================#
   """------------------------------------------------------------------
      init_struct_param: initialize structural parameters
   """
   def init_struct_param(self):
      self.nNodes = 0      
      self.nAtoms = 0
      self.nHterm = 0
      self.nRings = 0
      self.nAhold = 0
      self.nLocal = 0      
      self.nDeloc = 0

      self.lstAtSymb    = []
      self.lstAtNumb    = []
      self.lstAtLabel   = []
      self.lstAtCharge  = []
      self.lstAtDepth   = []     
      self.lstAromat    = []
      self.lstAaaEntr   = []      
      self.lstCaaEntr   = []
      self.lstAtf0      = []
      self.lstHcount    = []
      self.lstNbors     = []      
      self.lstNvbors    = []     
      self.lstRingMap   = []

      self.dictNbors = {}
      self.dictBonds = {}
      self.dictPairs = {}

      self.lstRingLength = []
      self.lstRingAromat = []
      self.lstRingCharge = []
      self.lstRings      = []

      self.lstIdxMat  = []
      self.lstAdjMat  = []
      self.lstDistMat = []

   #===================================================================#
   # RESET                                                             #
   #===================================================================#
   """------------------------------------------------------------------     
      reset: set new CurlySMILES notation and set parameters to
             initial (unassigned) values
             (a following call of parse() should assigns parameters for
              this new notation)
   """    
   def reset(sAnnSmi):
      self.sAnnSmi  = sAnnSmi
      self.lstTokens = None
      self.lstTokTyp = None
      self.lstDepth  = None
      self.maxDepth  = None       
      self.init_struct_param()
      self.lstErrors = []

   #===================================================================#
   # PREPARSE                                                          #
   #===================================================================#
   """------------------------------------------------------------------     
       preparse: preparse self.sAnnSmi and

       assignments: self.lstTokens, self.lstTokTyp, self.lstDepth,
                    maxDepth  

          TokTyp:  a = aromatic atomic symbol, b = bond symbol,
                   c = curly annotation, n = non-aromatic symbol,
                   q = square bracket notation, r = ring closure
                   
       EXAMPLE:
          self.sAnnSmi = 'o1c(Br)ccc1' 
          self.lstTokens = ['o','1','c','(','Br',')','c','c','c','1']
          self.lstTokTyp = ['a','r','a','(','n', ')','a','a','a','r']
          self.lstDepth  = [ 0,  0,  0,  1,  1,   1,  0,  0,  0,  0 ]
          self.maxDepth  = 1 

       return: 1 if successful, 0 otherwise with message in lstErrors
   """
   def preparse(self):

      self.lstTokens = []
      self.lstTokTyp = []
      self.lstDepth  = []
      self.maxDepth  = 0

      cTokTyp = '?'
      nDepth  = 0     # nesting depth
      nCurlyDepth = 0 # curly depth
      bInsideSqb  = 0 # inside square brackets
      sToken  = ''
      iChar = 0
      bSkip = 0
      for cChar in self.sAnnSmi:
         iChar += 1

         if bSkip != 0:
            bSkip = 0
            continue
         
         cLookAhead = '0'         
         if iChar < len(self.sAnnSmi):
            cLookAhead = self.sAnnSmi[iChar]
         
         if cTokTyp != '?':
            if (cTokTyp == 'c' and cChar == '}' and nCurlyDepth == 1)\
               or (cTokTyp == 'q' and cChar == ']'):
               self.lstTokens.append(sToken)
               self.lstTokTyp.append(cTokTyp)
               self.lstDepth.append(nDepth)      
               cTokTyp = '?'
               sToken  = ''
               nCurlyDepth = 0
               if cChar == ']':
                  bInsideSqb = 0
               continue # with next cChar
            elif cTokTyp == 'r':
               if cChar == '%' or cChar.isdigit():               
                  sToken += cChar
                  continue # with next cChar
               else:                    # done with ring closure(s)
                  self.lstTokens.append(sToken)
                  self.lstTokTyp.append(cTokTyp)
                  self.lstDepth.append(nDepth)      
                  cTokTyp = '?'
                  sToken  = ''
                  # continue with current cChar
            else:
               sToken += cChar
               if cChar == '{':
                  nCurlyDepth += 1
               elif cChar == '}':
                  nCurlyDepth -= 1                 
               continue # with next cChar

         if cChar == '{':
            cTokTyp = 'c'
            nCurlyDepth += 1
         elif cChar == '[':
            cTokTyp = 'q'
            bInsideSqb = 1
         elif cChar == '%' or cChar.isdigit():
            cTokTyp = 'r'
            sToken = cChar
         elif cChar == '(':
            nDepth += 1
            if nDepth > self.maxDepth:
               self.maxDepth = nDepth
            self.lstTokens.append(cChar)
            self.lstTokTyp.append(cChar)
            self.lstDepth.append(nDepth)
         elif cChar == ')':           
            self.lstTokens.append(cChar)
            self.lstTokTyp.append(cChar)
            self.lstDepth.append(nDepth)
            nDepth -= 1           
         elif cChar == '=' or cChar == '#' or cChar == '$':
            self.lstTokens.append(cChar)
            self.lstTokTyp.append('b')
            self.lstDepth.append(nDepth)            
         elif cChar == '*' or  cChar == '~':
            self.lstTokens.append(cChar)
            self.lstTokTyp.append(cChar)
            self.lstDepth.append(nDepth)                      
         elif cChar.isalpha():
            
            # check one-letter lower-case symbol (aromatic symbol)
            if cChar in ['b','c','n','o','s']:
               self.lstTokens.append(cChar) # one-letter, aromatic
               self.lstTokTyp.append('a')
               self.lstDepth.append(nDepth)
               continue
            elif not cChar.isupper():
               sMsg  = 'csm_subnotation.preparse: unxepected' 
               sMsg += ' lower-case letter near idx="%d"' % iChar
               self.lstErrors.append(sMsg)                           
               return 0

            # check upper-case symbol   
            if cLookAhead.isalpha() and cLookAhead.islower():
               sSymb2l = cChar + cLookAhead  # two-letter symbol ?
               if self.oDataFace.is_valid_symbol(sSymb2l) and \
                      self.is_2ch_atsymb(bInsideSqb,sSymb2l):
                  self.lstTokens.append(sSymb2l)
                  bSkip = 1 # skip next cChar in loop
               else:
                  self.lstTokens.append(cChar)
               self.lstTokTyp.append('n')
               self.lstDepth.append(nDepth)               
            elif self.oDataFace.is_valid_symbol(cChar): 
               self.lstTokens.append(cChar)
               self.lstTokTyp.append('n')
               self.lstDepth.append(nDepth)                       
            else:
               sMsg  = "preparse: invalid atomic " 
               sMsg += "symbol near iChar=%d" % iChar
               self.lstErrors.append(sMsg)               
               return 0

      # check if undone ring closure
      if cTokTyp == 'r':              
         self.lstTokens.append(sToken)
         self.lstTokTyp.append(cTokTyp)
         self.lstDepth.append(nDepth)      
         cTokTyp = '?'

      # final check on success, that is, consistency
      if cTokTyp != '?':
         sMsg  = "preparse: cTokTyp='%s'" % cTokTyp
         sMsg += " instead of '?' (unbalanced brackets?)"
         self.lstErrors.append(sMsg)
         return 0      
      elif nDepth != 0:
         sMsg  = "preparse: nDepth='%s'" % cTokTyp
         sMsg += " instead of '0' (unbalanced parenthesis?)"
         self.lstErrors.append(sMsg)
         return 0
      else: 
         return 1        

   """
       is_2ch_atsymb: check if a 2-char atomic symbol, which is not
                      enclosed in square brackets, needs to or
       otherwise represents two symbols, while the second lower-case
       symbol is an aromatic-ring atom.
       For example, notation CSc1ccccc1 (methyl pheny sulfide), does not
       contain scandium (Sc), Sc encodes S followed by aromatic C atom
       (Sc is not in the organic set and always has to be inside square
        brackets)                          
   """
   def is_2ch_atsymb(self, bInsideSqb, sSymb):
      if bInsideSqb == 0 and \
         (not self.oDataFace.is_in_organic_set(sSymb)):
         return 0
      else:
         return 1   

   #===================================================================#
   # EVALUATE: parse and generate molecular structure data             #
   #===================================================================#
   """------------------------------------------------------------------     
      parse: parse, assign member variables, including call of methods
             to generate distance matrix, atom fragment partition and
             molecular formula
      return: self.lstErrors
   """
   def parse(self):

      if self.lstTokens == None:
         self.preparse()
      self.init_struct_param()
      
      iToken = 0
      idxCurly = 0
      lstBrnchPnts  = [] # list of indexes of current branch point atoms,
                         # nDepth is index in lstBranchPoint
      dictOpenRings = {} # keeps track if open rings:
                         # key = ring id, value = [idx,lstBrnchPnts]
                         #                         of waiting atom
      dictRingData  = {} # {sAtPair: [lstBrnchPnts1,lstBrnchPnts2],...}
                         #  sAtPair = 'i,j' with i and j indexes of
                         #            ring closure atoms, i < j
                         #  lstBrnchPnts1 = branch points for atom i
                         #  lstBrnchPnts2 = branch points for atom j
      sWaitingBond = None                    
      for sToken in self.lstTokens:
         sTokTyp = self.lstTokTyp[iToken]
         nDepth  = self.lstDepth[iToken]    

         if sTokTyp == '(':
            pass
         elif sTokTyp == ')':   
            del lstBrnchPnts[len(lstBrnchPnts)-1]
         elif sTokTyp in ['b','~']: 
            sWaitingBond = sToken   
         elif sTokTyp == 'c': # annotation in curly braces
            oCurly = csm_curlyann.CurlyAnnotation(self.oDataFace,sToken)
            lstErrorsCurly = oCurly.parse()
            if len(lstErrorsCurly) > 0:
               for sMfErr in lstErrorsCurly:
                  sMsg = 'AnnotatedSmiles.parse: oCurly.%s' % sMfErr
                  self.lstErrors.append(sMsg)
            sAM = oCurly.annotation_marker()
            if self.oDataFace.is_aaa(sAM):  # if atom-anchored annotation
               (self.lstAaaEntr[self.nAtoms-1]).append(oCurly)              
            elif self.oDataFace.is_caa(sAM):  # if component-anchored annotation
               self.lstCaaEntr.append(oCurly)   
            else:
               sMsg = "AnnotatedSmiles.parse: unknown descriptor/annotation marker '%s'" % sAM
               self.lstErrors.append(sMsg)               
         elif sTokTyp == 'r':
            (lstRid,sErr) = self.parse_ring_closure(sToken)
            if sErr != None:
               sMsg  = 'AnnotatedSmiles.parse: %s' % sErr
               self.lstErrors.append(sMsg)             
               return  self.lstErrors
            else:
               for nRid in lstRid:
                  idxAt2 = self.nAtoms - 1                 
                  if dictOpenRings.has_key(nRid):
                     lstRingData = dictOpenRings[nRid]
                     idxAt1 = lstRingData[0]
                     sBond = '-'
                     if self.lstAromat[idxAt1] and self.lstAromat[idxAt2]:
                        sBond = ':'
                     del dictOpenRings[nRid]       
                     self.connect_new_atom(idxAt1,idxAt2,sBond)
                     # keep ring data for later evaluation
                     sAtPair = '%d,%d' % (idxAt1,idxAt2) 
                     lstBrnchPnts1 = lstRingData[1]
                     dictRingData[sAtPair] = [lstBrnchPnts1,lstBrnchPnts[:]]   
                  else:
                     dictOpenRings[nRid] = [idxAt2,lstBrnchPnts[:]]               
         elif sTokTyp in ['a','n','q','*']:
            self.nNodes += 1
            self.lstAaaEntr.append([])
            self.lstRingMap.append([])            
            idxCurly = 0            
            if sTokTyp == 'n' or sTokTyp == '*':
               self.lstAtSymb.append(sToken)
               nAtNumb = self.oDataFace.atnumb_as_int(sToken)
               if sToken == '*':
                  self.Ahold += 1
               if nAtNumb == None:
                  sMsg  = 'AnnotatedSmiles.parse: %d. atomic node' %\
                          (self.nAtoms+1)
                  sMsg += ' "%s",without atomic number' % sToken
                  self.lstErrors.append(sMsg)             
                  return  self.lstErrors
               self.lstAtNumb.append(nAtNumb)
               self.lstAtLabel.append(None)
               self.lstAtCharge.append(0)
               self.lstAtDepth.append(self.lstDepth[iToken])     
               self.lstAromat.append(0)
               self.lstHcount.append(None)  
            elif sTokTyp == 'a':
               sUpper = sToken.upper() 
               self.lstAtSymb.append(sUpper)
               nAtNumb = self.oDataFace.atnumb_as_int(sUpper)
               if nAtNumb == None:
                  sMsg  = 'AnnotatedSmiles.parse: %d.atomic' % self.nAtoms
                  sMsg += ' symbol, "%s",without atomic number' % sToken
                  self.lstErrors.append(sMsg)             
                  return  self.lstErrors
               self.lstAtNumb.append(nAtNumb)
               self.lstAtLabel.append(None)
               self.lstAtCharge.append(0)
               self.lstAtDepth.append(self.lstDepth[iToken])     
               self.lstAromat.append(1)
               self.lstHcount.append(None)                
            else:
               (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr) = \
                        self.parse_square_bracket_notation(sToken)               
               if sErr != None:
                  sMsg  = 'AnnotatedSmiles.parse: %d. atomic node' %\
                          (self.nAtoms + 1)
                  sMsg += ' < %s' % sErr
                  self.lstErrors.append(sMsg)             
                  return  self.lstErrors
               self.lstAtSymb.append(sAtSymb)
               nAtNumb = self.oDataFace.atnumb_as_int(sAtSymb)               
               self.lstAtNumb.append(nAtNumb)
               self.lstAtLabel.append(nLabel)
               self.lstAtCharge.append(int(sCharge))
               self.lstAtDepth.append(self.lstDepth[iToken])     
               self.lstAromat.append(bArom)
               self.lstHcount.append(nHAt)

            # connect current atom to left-side neighbor atom, idxNbor
            if self.nAtoms < 1:
               self.dictNbors[0] = []
               self.dictBonds[0] = []
            else:   
               sBond = '-'
               idxNbor = lstBrnchPnts[len(lstBrnchPnts)-1]
               if sWaitingBond != None:
                  sBond = sWaitingBond
                  sWaitingBond = None
               elif self.lstAromat[self.nAtoms] and \
                    self.lstAromat[idxNbor]:
                  sBond = ':'
               self.connect_new_atom(self.nAtoms,idxNbor,sBond)
            
            # update branch points and atom count   
            if len(lstBrnchPnts) <= nDepth:
               lstBrnchPnts.append(self.nAtoms) 
            else:
               lstBrnchPnts[nDepth] = self.nAtoms
            self.nAtoms += 1
         else:
            sMsg  = 'AnnotatedSmiles.parse: unknown token type ' 
            sMsg += '"%s" for token "%s",' % (sTokTyp,sToken)
            self.lstErrors.append(sMsg)             
            return self.lstErrors         

         iToken += 1

      # apply valence rule and derive count of non-node H-atoms
      self.count_tonode_valence_bonds()
      self.count_non_node_hatoms()
      self.nAtoms += self.nHterm

      # make adjacency matrix
      self.make_adj_mat()

      # make distance matrix
      self.make_dist_mat()

      # derive rings with member atoms
      self.make_rings(dictRingData)

      # shorten rings
      self.shorten_rings()

      # revisit rings to assign remaining ring data
      self.revisit_rings()

      # sum charges
      self.sum_local_charges()      
      self.sum_deloc_charges()

      # make molecular formula
      self.make_mf()

      # make atom fragments
      self.make_atf0()

      return self.lstErrors

   """------------------------------------------------------------------   
      parseSquareBracketNotation: parse sStr with notation from square
                                  bracket enclosure

      return: (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)
               nLabel  = isotope label (positive integer)
               sAtSymb = atomic symbol, capitalized if encountered as
                         lower-case aromatic symbol
               bArom   = 0 = non-aromatic, 1 = aromatic ring atom          
               nHat    = number of adjacent H-atoms, 0 if none
               sCharge = '0' for neutral, 'n' for positive, '-n' for
                         negative charge
               sErr    = None, if successful, error message otherwise          
              
      Example 1: sStr = '15NH4+'
                 returns: (15,'N',0, 4, '1', None)
               
      Example 2: sStr = 'si'  # Si-atom in aromatic ring
                 returns: (None, 'Si', 1, 0, '0', None)
               
   """
   def parse_square_bracket_notation(self,sStr):

      nLabel  = None
      sAtSymb = None
      bArom   = None
      nHAt    = None
      sCharge = None
      sErr    = None

      sRemaining = sStr

      # check for positive charge value
      lstPair = sRemaining.split('+',1)
      if len(lstPair) == 2:
         sFound =  lstPair[1]
         sRemaining = lstPair[0]
         if len(sFound) == 0:
            sCharge = '1'  # +1
         elif sFound.isdigit():
            sCharge = sFound
         else:
            # check for multiple + signs
            cntPlus = sFound.count('+')
            if cntPlus == len(sFound):
               sCharge = '%d' % (cntPlus + 1)
            else:
               sErr  = "parse_square_bracket_notation(%s):" %\
                        sStr         
               sErr += 'with miscoded positive charge'  
               return (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)
                  
      # check for negative charge value, if no positive was found
      if sCharge == None:
         lstPair = sRemaining.split('-',1)
         if len(lstPair) == 2:
            sFound =  lstPair[1]
            sRemaining = lstPair[0]
            if len(sFound) == 0:
               sCharge = '-1' 
            elif sFound.isdigit():
               sCharge = '-%s' % sFound
            else:
               # check for multiple minus signs
               cntMinus = sFound.count('-')
               if cntMinus == len(sFound):
                  sCharge = '-%d' % (cntMinus + 1)
               else:
                  sErr  = "parse_square_bracket_notation(%s):" %\
                       sStr         
                  sErr += 'with miscoded negative charge'  
                  return (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)
         else:                  
            sCharge = '0'

      # check for isotope label
      sLabel = ''
      for cChar in sRemaining:
         if cChar.isdigit():
            sLabel += cChar
         else:
            break   
      lenLabel = len(sLabel)
      if lenLabel > 0:
         nLabel = int(sLabel)
         sRemaining = sRemaining[lenLabel:]
      
      # get atom symbol, the only required part inside []
      if len(sRemaining) < 1:
         sErr  = "parse_square_bracket_notation(%s):" %\
                sStr         
         sErr += 'missing atom symbol'  
         return (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)
      elif len(sRemaining) == 1:
         sAtSymb = sRemaining
      else:   
         sAtSymb = sRemaining[0:2]
      # verify atom symbol
      if sAtSymb[0] == '*':
         sAtSymb = '*'
      else:   
         if not sAtSymb.isalpha():
            sErr  = "parse_square_bracket_notation(%s):" %\
                   sStr         
            sErr += 'atomic symbol contains non-letter(s)'  
            return (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)
      # check for aromatic atom, capitalize first letter if aromatic
      cChar = sAtSymb[0]
      if cChar.islower():
         bArom = 1
         if len(sAtSymb) > 1:
            sAtSymb = cChar.upper() + sAtSymb[1]
         else:
            sAtSymb = cChar.upper()                
      else:
         bArom = 0
      # check if valid atom symbol   
      if not self.oDataFace.is_valid_symbol(sAtSymb):   
         if len(sAtSymb) > 1: # check
            sAtSymb = sAtSymb[0]
            if not self.oDataFace.is_valid_symbol(sAtSymb):            
              sErr  = "parse_square_bracket_notation(%s): " %\
                sStr         
              sErr += 'invalid atomic symbol'  
              return (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)

      if sAtSymb == 'H' and len(sRemaining) > 1:
         sErr  = "parse_square_bracket_notation(%s):" %\
                sStr         
         sErr += '"featured H-atom" cannot be followed by H-counting code'  
         return (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)
             
      # check for H-atom(s)
      if len(sRemaining) > len(sAtSymb):
         sRemaining = sRemaining[len(sAtSymb):]
         if sRemaining[0] != 'H':
            sErr  = "parse_square_bracket_notation(%s):" %\
                sStr         
            sErr += '"H" (or charge sign) expected to follow atomic symbol'  
            return (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)
         if len(sRemaining) == 1:
            nHAt = 1
         elif len(sRemaining) > 1:   
            if sRemaining[1:].isdigit():
               nHAt = int(sRemaining[1:])
            else:
               sErr  = "parse_square_bracket_notation(%s):" %\
                   sStr         
               sErr += 'syntax error after H'        
      else:
         nHAt = 0
      
      return (nLabel,sAtSymb,bArom,nHAt,sCharge,sErr)

   """------------------------------------------------------------------ 
      parse_ring_closure: parse sStr with notation for the closure
                        of one or more rings (sStr is expected to
                        contain only digits and %)

      return: (lstRid,sErr)
               lstRid  = list of ring identifiers (whole integers)
               sErr    = None, if successful, error message otherwise  
              
      Example 1: sStr = '21'
                 returns: ([2,1],None)
               
      Example 2: sStr = '56%12'  # 
                 returns: ([5,6,12],None)
               
   """
   def parse_ring_closure(self,sStr):

      lstRid  = [] # list of ring identifiers 
      sErr    = None

      if sStr.isdigit():
         for cDigit in sStr:
            lstRid.append(int(cDigit))
      else:
         lstParts = sStr.split('%')
         iPart = 0
         for sPart in lstParts:
            iPart += 1             
            if iPart == 1:
               if len(sPart) < 1:
                  continue
               for cChar in sPart:
                 if cChar.isdigit():
                    lstRid.append(int(cChar))
                 else:
                    sErr  = "parse_ring_closure(%s):" %\
                    sStr         
                    sErr += 'unexpected char "%s"' % cChar
                    return (lstRid,sErr) 
            else:
               if sPart.isdigit():
                  lstRid.append(int(sPart))
               else:
                  sErr  = "parse_ring_closure(%s):" %\
                  sStr         
                  sErr += 'invalid ring identifier: "%s"' % sPart
                  return (lstRid,sErr)                
             
      return (lstRid,sErr)

   #===================================================================#
   # CONNECTIVITY                                                      #
   #===================================================================#
   """------------------------------------------------------------------    
      connect_new_atom: connect a new atom (idxAt1) with a previously
                        entered atom (idxAt2) by bond given with sBond
      update: self.dictNbors, self.dictBons                   
   """
   def connect_new_atom(self,idxAt1,idxAt2,sBond):

      if self.dictNbors.has_key(idxAt1):
         lst = self.dictNbors[idxAt1]
         lst.append(idxAt2)
         self.dictNbors[idxAt1] = lst
         lst = self.dictBonds[idxAt1]
         lst.append(sBond)
         self.dictBonds[idxAt1] = lst
      else:   
         self.dictNbors[idxAt1] = [idxAt2]      
         self.dictBonds[idxAt1] = [sBond]  

      if self.dictNbors.has_key(idxAt2):
         lst = self.dictNbors[idxAt2]
         lst.append(idxAt1)
         self.dictNbors[idxAt2] = lst
         lst = self.dictBonds[idxAt2]
         lst.append(sBond)
         self.dictBonds[idxAt2] = lst
      else:   
         self.dictNbors[idxAt2] = [idxAt1]      
         self.dictBonds[idxAt2] = [sBond]  

      sPair = '%d,%d' % (idxAt1,idxAt2)
      self.dictPairs[sPair] = sBond
      sPair = '%d,%d' % (idxAt2,idxAt1)
      self.dictPairs[sPair] = sBond


   """------------------------------------------------------------------    
      count_tonode_valence_bonds: count the number of valence bonds
                                  that connect a node atom to its
                                  neigbor node atoms
                                  
      EXAMPLES: 1. chloroethene: ClC=C
                                 idxAt = 0: Cl has 1 tonode valence bond
                                 idxAt = 1: C has 3 tonode valence bons 
                                 idxAt = 2: C has 2 tonode valence bonds

                2. pyridine: c1nc1ccn1                            
                             each C or N atom has 3 valence bonds
      
      assign: self.lstNbors and self.lstNvbors         
   """
   def count_tonode_valence_bonds(self):

      """
          dictionary with bond-symbol/valence-count values
          (values are doubled to avoid floating point numbers for
           bonds in aromatic rings)  
      """
      dictVB = {'-' : 2, '=' : 4, '#' : 6, '$' : 8, ':' : 3}

      for idxAt in range(self.nAtoms):
         lstBonds = self.dictBonds[idxAt] 
         lstAtCrl = self.lstAaaEntr[idxAt]

         self.lstNbors.append(len(lstBonds))  # number of neigbor nodes

         nVbors = 0 # valence bonds connecting to neighbors
         for sBond in lstBonds: # for each bond to a neighbor node
            
            if dictVB.has_key(sBond):
               nVbors += dictVB[sBond]
            else:
               nVbors = None
               break
            
         # count open valence bonds given in annotations
         if nVbors != None:
            for oCrl in lstAtCrl:
               sAM  = oCrl.annotation_marker()
               nCor = self.oDataFace.hcount_correction(sAM)
               nVbors += nCor                        

         # update self.lstNvbors
         if nVbors != None:      
            nVbors /= 2
            self.lstNvbors.append(nVbors)
         else:
            self.lstNvbors.append(0)


   """------------------------------------------------------------------    
      count_non_node_hatoms: count the non-node hydrogen atoms that
                             are adjacent to a node atom

      apply_valence_rule: apply valence rule, i.e.  normal valence
                          assumption:
         (a) single normal valence of B, C, O, and halogens
             are                      3, 4, 2, and 1;
         (b) for nitrogen and phosphorus 3 or 5;
         (c) for aliphatic sulfur 2, 4, or 6;
         (d) atoms in other states and  non-organic-subset atoms
             have explicit H-counts given in square brackets

      assign: self.lstHcount and self.nHterm
   """
   def count_non_node_hatoms(self):
      
      for idxAt in range(self.nAtoms):
         nAtNumb  = self.lstAtNumb[idxAt]
         sAtSymb  = self.lstAtSymb[idxAt]
         nNbors   = self.lstNbors[idxAt]
         nNvbors  = self.lstNvbors[idxAt]

         if self.lstHcount[idxAt] != None:
            self.nHterm += self.lstHcount[idxAt]
         else:   
            nDiff = 0 # difference between expected and obtained
                      # values for number of valence bonds
            if nAtNumb == 5:               # boron
               nDiff = 3 - nNvbors
            elif nAtNumb == 6:             # carbon
               nDiff = 4 - nNvbors              
            elif nAtNumb == 9:             # fluorine
               nDiff = 1 - nNvbors
            elif nAtNumb == 17:            # chlorine
               nDiff = 1 - nNvbors
            elif nAtNumb == 35:            # bromine
               nDiff = 1 - nNvbors
            elif nAtNumb == 53:            # iodine
               nDiff = 1 - nNvbors
            elif nAtNumb == 7:             # nitrogen
               if self.lstAromat[idxAt]:
                  if nNvbors >= 3:   # such as in pyridine-N-oxide 
                     nDiff = 0       # or in N-substituted 5-ring
                  else:
                     nDiff = 0       # in 6-ring
               else:   
                  if nNvbors > 3:
                     nDiff = 5 - nNvbors
                  else:
                     nDiff = 3 - nNvbors
            elif nAtNumb == 8:             # oxygen
               if self.lstAromat[idxAt]:
                  nDiff = 0
               else:
                  nDiff = 2 - nNvbors              
            elif nAtNumb == 15:            # phosphorus
               if self.lstAromat[idxAt]:
                  if nNvbors >= 4:
                     nDiff = 0       # in P-substituted 5-ring
                  else:
                     nDiff = 0       # in 6-ring
               elif nNbors == 0:
                     nDiff = 3       # PH3
               elif nNbors == 1:
                  lstNbors = self.dictNbors[idxAt]
                  nAtNumbNbor = self.lstAtNumb[lstNbors[0]] 
                  if nNvbors == 1:
                     nDiff = 2       # H2P-PH2, R-PH2, X-PH2
                  elif nNvbors == 2:
                     if nAtNumbNbor in [8,16]:
                        nDiff = 3    # O=PH3, S=PH3 ?
                     else:
                        nDiff = 1    # >C=PH ?
               elif nNbors == 2:
                  if nNvbors == 2:
                     nDiff = 1       # H2P-PH-PH2; cyclic P3H3, P4H4; R-PH-R
                  elif nNvbors == 3:
                     nDiff = 2       # O=PH2OH hypophosphorous acid
                                     #         (phosphinic acid) and derivatives
               elif nNbors == 3:
                  if nNvbors == 4:
                     nDiff = 1       # O=PH(OH)2 and derivatives
               # else  any neutral PH, PH2  with 3 or more neighbors?  

            elif nAtNumb == 16:             # sulfur
               if self.lstAromat[idxAt]:
                  nDiff = 0
               elif nNbors == 0:
                  nDiff = 2                 # H2S
               elif nNbors == 1:
                  lstNbors = self.dictNbors[idxAt]
                  nAtNumbNbor = self.lstAtNumb[lstNbors[0]] 
                  if nNvbors == 1:
                     nDiff = 1              # example: HS-R 
                  # elif nNvbors == 2:
                  #   if nAtNumbNbor in [6, 7, 14]:
                  #      nDiff = 0           # example: S=C<, S=N-, S=P-,
                  #   else:
                  #      nDiff = 0           # example: H2S=O                      
               elif nNbors == 2:
                  if nNvbors == 2:          # example: -S-
                     nDiff = 0
                  elif nNvbors == 3:
                     nDiff = 1              # example: -SH=
               elif nNbors == 3:
                  if nNvbors == 5:          # example  O=SH(CH3)=O
                     nDiff = 1                  
               elif nNbors == 4:
                  if nNvbors < 6:
                     nDiff = 6 - nNvbors    
            else:
               sMsg  = 'csm_subnotation.applyValenceRule: "%s"' % sAtSymb
               sMsg += ' not in organic subset, hence explicit Hcount'
               sMsg += ' assignment expected' 
               self.lstErrors.append(sMsg)             
               return 0            
            
            self.lstHcount[idxAt] = nDiff
            self.nHterm += nDiff

   #===================================================================#
   # TOPOLOGICAL MATRICES                                              #
   #===================================================================#
   """------------------------------------------------------------------       
      make_adj_mat: make adjacency matrix:

                  lstIdxMat = [(0, 1), (0, 2),...,(nNodes-2,nNodes-1)]
                  lstIdxMat is a list of pairs (idxRow,idxCol) to
                            look up  
      k = 0 # index for lstIdxMat, k=0,1,2,...,nNodes*(nNodes-1)/2                   
   """
   def make_adj_mat(self):

      self.lstIdxMat  = []
      nRowLength = self.nNodes-1  # length of row in upper triangle
      for idxRow in range(self.nNodes-1):
         for idxCol in range(nRowLength): 
            self.lstIdxMat.append( (idxRow,idxCol+idxRow+1) )
         nRowLength -= 1   

      self.lstAdjMat = []
      for (i,j) in self.lstIdxMat:
         if self.are_connected(i,j):
            self.lstAdjMat.append(1)
         else:   
            self.lstAdjMat.append(0)
      
   """------------------------------------------------------------------                
      make_dist_mat: make topological distance matrix:
                   based on lstAdjMat, calculate lstDistMat using
                   a published algorithm:

                   W. R. Mueller, K. Szymanski, J. V. Knop and
                   N. Trinajstic:
                   An Algorithm for Construction of the Molecular
                   Distance Matrix.
                   J. Comp. Chem. 1987, 8(2), pp. 170-173.
                   Table I therein contains the algorithm in FORTRAN.

      Example: 2-methylbutane, entered as CC(C)CC:
      
                   adjacency       distance 
                   matrix          matrix

                 i:0 1 2 3 4     i:0 1 2 3 4              
               j:    
               0   0 1 0 0 0       0 1 2 2 3
               1   1 0 1 1 0       1 0 1 1 2
               2   0 1 0 0 0       2 1 0 2 3
               3   0 1 0 0 1       2 1 2 0 1
               4   0 0 0 1 0       3 2 3 1 0
               
               The serialized upper triangle of each matrix:
               lstAdjMat  = [1,0,0,0,1,1,0,0,0,1]
               lstDistMat = [1,2,2,3,1,1,2,2,3,1]
               idxK(i,j):    0 1 2 3 4 5 6 7 8 9                   
   """
   def make_dist_mat(self):

      # no distance matrix if no adjacency matrix
      if self.lstAdjMat == None: 
         return None         

      # initialize distance matrix
      self.lstDistMat = []
      for nEntry in self.lstAdjMat:
         if nEntry == 1:   
            self.lstDistMat.append(1)
         else:
            self.lstDistMat.append(self.nAtoms)

      # algorithm of Mueller, Szymanski, Knop, and Trinajstic
      L = 1
      while L < self.nNodes-1:
         for J in range(self.nNodes):
            for I in range(self.nNodes):
               if I==J:
                  continue
               idxK = self.idxK(I,J)
               if self.lstDistMat[idxK] == L:
                  for K in range(self.nNodes):
                     idxKI = self.idxK(K,I)
                     iakj = self.lstDistMat[idxKI] + L
                     idxKJ = self.idxK(K,J)
                     if self.lstDistMat[idxKJ] > iakj:
                        self.lstDistMat[idxKJ] = iakj
         L += L
         
      return self.lstDistMat

   """------------------------------------------------------------------ 
      idxK: return index k of serialized upper triangle for entry i,j
            in square symmetric adjacency or distance matrix 
   """
   def idxK(self,i,j):
      if i >= self.nNodes or j >= self.nNodes or \
         i < 0 or j < 0:
         sErr  = "csm_subnotation.idxK(%d,%d): exceeding " % (i,j)  
         sErr += 'allowed index range for nAtoms=%d' % self.nNodes
         self.lstErrors.append(sErr)
         return None
      elif i==j:
         return -1
      elif i < j:
         ii = j
         jj = i
      else:
         ii = i
         jj = j
      k = 0
      for kk in range(jj):
         k += (self.nNodes-1) - kk
      return (k+ii-jj)-1

   """------------------------------------------------------------------ 
      return topological distance between atoms i and j for any values
             of i and j between 0 and nAtoms-1
   """
   def dist(self,i,j):
      if i==j:
         return 0
      else:
         k = self.idxK(i,j)
         if k != None:         
            return self.lstDistMat[k]
         else:
            return None

   #===================================================================#
   # EVALUATE RINGS                                                    #
   #===================================================================#
   """------------------------------------------------------------------    
      make_rings: recognize rings based on dictRingData with format:

      dictRingData  = sAtPair: [lstBrnchPnts1,lstBrnchPnts2],...}
                         sAtPair = 'i,j' with i and j indexes of
                                     ring closure atoms, i < j
                         lstBrnchPnts1 = branch points for atom i
                         lstBrnchPnts2 = branch points for atom j
      
      assign: self.lstRings
   """
   def make_rings(self,dictRingData):

      self.lstRings = []

      for sAtPair in dictRingData.keys():
         lstPair = sAtPair.split(',')
         idxAt1 = int(lstPair[0])
         idxAt2 = int(lstPair[1])
         [lstBrnchPnts1,lstBrnchPnts2] = dictRingData[sAtPair]           
         nDepth1 = self.lstAtDepth[idxAt1] 
         nDepth2 = self.lstAtDepth[idxAt2] 

         lstMembers = []
         if nDepth1 == 0 and nDepth2 == 0:
            for k in range((idxAt2-idxAt1)+1):
                idxBetween = idxAt1 + k
                if self.lstAtDepth[idxBetween] == 0:
                   lstMembers.append(idxBetween) 
         else:
            # minDepth = min(nDepth1,nDepth2)
            # For each ring closing get ring members down to branch point
            # with depth value equal to minDepth
            [idxLeftMost,lstRight] = \
            self.ring_members_right_to_left(idxAt1,idxAt2,lstBrnchPnts1,lstBrnchPnts2)
#            print 'idxAt1=',idxAt1,'lstBrnchPnts1:', lstBrnchPnts1
#            print 'idxAt2=',idxAt2,'lstBrnchPnts2:', lstBrnchPnts2            
#            print '### lstRight:', lstRight
            lstRight.sort()
#            print '#### idxLeftMost=',idxLeftMost,', idxAt1=',idxAt1
            
            if idxLeftMost == idxAt1:
               lstMembers = lstRight               
            elif idxLeftMost in lstBrnchPnts1:
               lstLeft = self.ring_members_left_to_right(idxLeftMost,idxAt1)            
               lstMembers = lstLeft + lstRight
            else:
               lstMembers = [idxAt1] + lstRight

         self.lstRings.append(lstMembers)
         
   """------------------------------------------------------------------
      ring_members_right_to_left:
      return: [idxLeftMost,lstRight]
               idxLeftMost = index of left most atom that is ring rember
               lstRight    = list of ring members between idxLeftMost and idxAt2
   """
   def ring_members_right_to_left(self,idxAt1,idxAt2,lstBrnchPnts1,lstBrnchPnts2):
      idxLeftMost = -1
      lstRight = [idxAt2]
      nDepthMemb = len(lstBrnchPnts2)-1

      idxPrevMemb = idxAt2 # previous ring member atom
      for k in range(1,idxAt2+1):  
         idxAt = idxAt2 - k
         if self.are_connected(idxAt,idxPrevMemb):
            lstRight.append(idxAt)
            if idxAt == idxAt1 or idxAt in lstBrnchPnts1:
               idxLeftMost = idxAt
               break
            idxPrevMemb = idxAt

      return [idxLeftMost,lstRight]

   """------------------------------------------------------------------
      ring_members_left_to_right
   """
   def ring_members_left_to_right(self,idxLeftMost,idxAt1):

      lstLeft = [idxAt1]

      idxPrevMemb = idxAt1 # previous ring member atom
      for k in range(1,idxAt1):
         idxAt = idxAt1 - k
         if idxAt == idxLeftMost:
            break
         elif self.are_connected(idxAt,idxPrevMemb):
            lstLeft.append(idxAt)
         idxPrevMemb = idxAt

      return lstLeft

   """------------------------------------------------------------------ 
      shorten_rings: search for short-cuts in current rings
      update: self.lstRings
   """
   def shorten_rings(self):

      idxRing = 0
      for lstMemb in self.lstRings:
         lenRing = len(lstMemb)
         if lenRing > 3:

            lstMembCur = lstMemb[:]
            for k in range(lenRing-1):
               lstMembNew = self.get_next_shortcut(lstMembCur)
               if len(lstMembNew) < len(lstMembCur):
                  self.lstRings[idxRing] = lstMembNew
               else:
                  break

               lstMembCur = lstMembNew[:]    

         idxRing += 1
                  
   """------------------------------------------------------------------
      get_next_shortcut
      return: return list with ring members of shortcut ring
   """
   def get_next_shortcut(self,lstMemb):
      idxLast = len(lstMemb) - 1  # index of last atom in sequence

      idxMemb1 = 0
      for idxAt1 in lstMemb[idxMemb1:idxLast-1]:

         idxMemb2 = idxMemb1+2
         for idxAt2 in lstMemb[idxMemb1+2:]:

            nDistInRing  = idxMemb2 - idxMemb1
            nDistFromMat = self.dist(idxAt1,idxAt2)

            if nDistFromMat < nDistInRing:

               lstShortcut = []
               lstShortPaths = self.shortest_paths(idxAt1,idxAt2)
               if len(lstShortPaths) < 1:
                  pass # ERROR
               else:
                  if nDistFromMat > 1:
                     lstShortcut = lstShortPaths[0]

               cntNew = 0
               for idxAtTrial in lstShortcut:
                  if not idxAtTrial in lstMemb:
                     cntNew += 1

               isHeadTailPair = 0
               if (idxAt1 == lstMemb[0] and idxAt2 == lstMemb[idxLast])  or \
                  (idxAt2 == lstMemb[0] and idxAt1 == lstMemb[idxLast]):
                  isHeadTailPair = 1

               if cntNew > 0:
                  lstSmallerRing = lstMemb[0:idxMemb1+1]
                  lstSmallerRing += lstShortcut
                  lstSmallerRing += lstMemb[idxMemb2:]   
                  return (lstSmallerRing)
               # elif idxAt1 and idxAt2 are adjacent, but not a head-tail pair
               elif  nDistFromMat == 1 and isHeadTailPair == 0:
                  lstSmallerRing = lstMemb[0:idxMemb1+1]
                  lstSmallerRing += lstMemb[idxMemb2:]   
                  return (lstSmallerRing)               
      
            idxMemb2 += 1

         idxMemb1 += 1
                    
      return (lstMemb[:])

   """------------------------------------------------------------------
      revisit_rings: over all node atoms, sum charge values 
      assign: self.lstRingLength, self.RingAromat, self.lstRingCharge   
   """   
   def revisit_rings(self):

      self.nRings = len(self.lstRings) 
      for iRing in range(self.nRings):

         lstMemb = self.lstRings[iRing]

         # ring length = number of ring member atoms
         self.lstRingLength.append(len(lstMemb))

         # flag aromatic rings and update atom-ring mapping
         isAromat = 1
         for iAtom in lstMemb:
            if self.lstAromat[iAtom] != 1: 
               isAromat = 0
            self.lstRingMap[iAtom].append(iRing)
            
         self.lstRingAromat.append(isAromat)            

         # ring charge: sum of local atom charges plus
         #              charges given with key e in annotations
         nLocal = self.local_ring_charge(iRing)
         nDeloc = self.deloc_ring_charge(iRing)
               
         self.lstRingCharge.append(nLocal+nDeloc)

   """------------------------------------------------------------------
      local_ring_charge: sum localized charges over all atom that are
                         a member of ring iRing
      return: nDeloc   
   """   
   def local_ring_charge(self, iRing):
      nLocal = 0
      lstMemb = self.lstRings[iRing]
      for iAtom in lstMemb:
         nLocal += self.lstAtCharge[iAtom]
      return nLocal

   """------------------------------------------------------------------
      deloc_ring_charge: calculate delocalized charge of ring iRing by
                         getting every e-value from annotation
      dictionaries that are anchored to a ring member atom
      return: nDeloc   
   """   
   def deloc_ring_charge(self, iRing):

      nDeloc = 0
      lstMemb = self.lstRings[iRing]
      lstAnchors = self.aaa_indices()  
      for iAtom in lstMemb:
         if iAtom in lstAnchors:
            lstCurlies = self.lstAaaEntr[iAtom]
            for oCurly in lstCurlies:          
               dictAnn = oCurly.annotation_dict()
               for sKey in dictAnn.keys():
                  if sKey == 'e':
                     sCharge = dictAnn[sKey]
                     if sCharge == '+':
                        nDeloc += 1
                     elif sCharge == '-':
                        nDeloc -= 1
                     elif len(sCharge) > 1:
                        if sCharge[0] == '+':
                           if sCharge[1:].isdigit:
                              nDeloc += int(sCharge[1:])
                           else:
                              sMsg  = 'deloc_ring_charges: atom %d ' % iAtom
                              sMsg += 'invalid positive charge: e=%s' % sCharge 
                              self.lstErrors.append(sMsg)            
                        elif sCharge[0] == '-':
                           if sCharge[1:].isdigit:
                              nDeloc -= int(sCharge[1:])
                           else:                         
                              sMsg  = 'deloc_ring_charges: atom %d ' % iAtom
                              sMsg += 'invalid negative charge: e=%s' % sCharge 
                              self.lstErrors.append(sMsg)          
                     else:
                        sMsg  = 'deloc_ring_charge: atom %d ' % iAtom
                        sMsg += 'invalid charge: e=%s' % sCharge 
                        self.lstErrors.append(sMsg)
                        
      return nDeloc

   #===================================================================#
   # CHARGE CONTRIBUTIONS                                              #
   #===================================================================#
   """
      sum_local_charges: over all node atoms, sum charge values 
      assign: self.nLocal   
   """   
   def sum_local_charges(self):
      for nCharge in self.lstAtCharge:
         self.nLocal += nCharge
   
   """
      sum_deloc_charges: over all node atoms, sum values given by
                         key e in annotation dictionaries
      assign: self.nDeloc   
   """   
   def sum_deloc_charges(self):

      iAtom = 0
      for lstCurlies in self.lstAaaEntr:
         for oCurly in lstCurlies:
            # sAM     = oCurly.annotation_marker()
            dictAnn = oCurly.annotation_dict()
            for sKey in dictAnn.keys():
               if sKey == 'e':
                  sCharge = dictAnn[sKey] 
                  if sCharge == '+':
                     self.nDeloc += 1
                  elif sCharge == '-':
                     self.nDeloc -= 1
                  elif len(sCharge) > 1:
                     if sCharge[0] == '+':
                        if sCharge[1:].isdigit:
                           self.nDeloc += int(sCharge[1:])
                        else:
                           sMsg  = 'sum_deloc_charges: atom %d ' % iAtom
                           sMsg += 'invalid positive charge: e=%s' % sCharge 
                           self.lstErrors.append(sMsg)            
                     elif sCharge[0] == '-':
                        if sCharge[1:].isdigit:
                           self.nDeloc -= int(sCharge[1:])
                        else:                         
                           sMsg  = 'sum_deloc_charges: atom %d ' % iAtom
                           sMsg += 'invalid negative charge: e=%s' % sCharge 
                           self.lstErrors.append(sMsg)          
                     else:
                        sMsg  = 'sum_deloc_charges: atom %d ' % iAtom
                        sMsg += 'invalid charge: e=%s' % sCharge 
                        self.lstErrors.append(sMsg)            
         iAtom += 1


   #===================================================================#
   # MOLECULAR FORMULA                                                 #
   #===================================================================#
   """------------------------------------------------------------------
      make_mf: generate molecular formula via dictOfDict
      assign:  self.oMf, an object of class csm_molform.MolecularFormula           
   """
   def make_mf(self):
     
      # non-node hydrogen atoms
      dictOfDict = { }
      if self.nHterm > 0:
         dictOfDict['H'] = { 'H': self.nHterm} 
              
      # 
      idxAt = 0
      for sAtSymb in self.lstAtSymb:
         nLabel   = self.lstAtLabel[idxAt]
         sLblSymb = sAtSymb
         if nLabel > 0:
            sLblSymb = '^%d%s' % (nLabel,sAtSymb)            

         if dictOfDict.has_key(sAtSymb):
            subdict = dictOfDict[sAtSymb]
            if subdict.has_key(sLblSymb):
               subdict[sLblSymb] = subdict[sLblSymb] + 1
            else:
               subdict[sLblSymb] = 1
         else:
            dictOfDict[sAtSymb] = {sLblSymb: 1}

         idxAt += 1

#      print dictOfDict

      nCharge = self.nLocal + self.nDeloc
      sCharge = '0'
      if nCharge > 0:
         sCharge = '%d+' % nCharge
      elif nCharge < 0:
         sCharge = '%d-' % ( (-1) * nCharge )
         
      self.oMf = csm_molform.MolecularFormula(self.oDataFace)
      self.oMf.set_dict_of_dict(dictOfDict,sCharge)
      lstMfErr = self.oMf.evaluate()
      if len(lstMfErr) > 0:
         for sMfErr in lstMfErr:
            sMsg = 'make_mf: oMf.%s' % sMfErr
            self.lstErrors.append(sMsg)

   #===================================================================#
   # ATOM FRAGMENTS                                                    #
   #===================================================================#
   """------------------------------------------------------------------
      make_atf0: generate atom fragments (level 0: atomic symbol of node
                 atom plus adjacent non-node H-atoms plus bonds to
                 neighbor nodes).
                 Each atf0 is a CurlySMILES notation.
      assign:  self.lstAtf0           
   """
   def make_atf0(self):

      idxAt = 0
      for sAtSymb in self.lstAtSymb:

         bSquareBrack = 0 # specifies if SQC encoding is required

         # get dangling bonds
         sSingleBond = ''
         sDoubleBond = ''
         sTripleBond = ''
         sQuapleBond = ''
         sAromatBond = ''
         sUnspecBond = ''
         lstBonds = self.dictBonds[idxAt] 
         for sBond in lstBonds:
            if sBond == '-':
               sSingleBond += '{-}'
            elif sBond == '=':
               sDoubleBond += '{=}'               
            elif sBond == '#':
               sTripleBond += '{#}'              
            elif sBond == '$':
               sQuapleBond += '{$}'
               bSquareBrack = 1
            elif sBond == ':':
               sAromatBond += '{:}'
               bSquareBrack = 1
            elif sBond == '~':
               sUnspecBond += '{~}'              
               bSquareBrack = 1

         # look for open bond contributions in annotations 
         lstCurlies  = self.lstAaaEntr[idxAt]
         for oCurly in lstCurlies:
            sAM = oCurly.annotation_marker()
            if self.oDataFace.is_open_single_bond(sAM):
               sSingleBond += '{-}'
            elif self.oDataFace.is_open_double_bond(sAM):
               sDoubleBond += '{=}'
            elif self.oDataFace.is_open_triple_bond(sAM):
               sTripleBond += '{#}'               
            elif self.oDataFace.is_open_quadruple_bond(sAM):
               sQuapleBond += '{$}'
               bSquareBrack = 1                
            elif self.oDataFace.is_open_aromatic_bond(sAM):
               sAromatBond += '{:}'
               bSquareBrack = 1                
            elif self.oDataFace.is_open_unspecified_bond(sAM):
               sUnspecBond += '{~}'
               bSquareBrack = 1                

         # get needed atom attributes
         nLabel  = self.lstAtLabel[idxAt]
         nCharge = self.lstAtCharge[idxAt]
         bAromat = self.lstAromat[idxAt] 
         nHcount = self.lstHcount[idxAt]

         # generate atomic symbol notation (include label,hcount,charge)
         if (not self.oDataFace.is_in_organic_set(sAtSymb)) or \
            len(sUnspecBond) > 0:
            bSquareBrack = 1
         
         if bAromat == 1:           
            sAtSymb = sAtSymb.lower()

         sAtf0 = sAtSymb
         if nLabel > 0:
             bSquareBrack = 1
             sAtf0 = '%d%s' % (nLabel,sAtf0)
         sCharge = ''
         if nCharge != 0:
             bSquareBrack = 1
             if nCharge == 1:
                sCharge = '+'
             elif nCharge == -1:
                sCharge = '-'
             elif nCharge > 1:
                sCharge = '+%d' % nCharge 
             elif nCharge < -1:   
               sCharge = str(nCharge) 

         if  bSquareBrack == 1:
            sAtf0 = '[' + sAtf0
            if nHcount > 0:
               sAtf0 += 'H'
               if nHcount > 1:
                  sAtf0 += str(nHcount)
            sAtf0 += sCharge + ']'   


         # add notation with open bond(s)
         sAtf0 += sAromatBond + sQuapleBond + sTripleBond + \
                  sDoubleBond + sSingleBond + sUnspecBond

         self.lstAtf0.append(sAtf0)

         idxAt += 1

   #===================================================================#
   # COMMON REQUESTS                                                   #
   #===================================================================#
   """------------------------------------------------------------------
      are_connected: are atoms with index idxAt1 and idxAt2 connected?
                    (idxAt1 and idxAt2: any int from 0 to self.nNodes-1)
      return: 1 or 0, depending on whether they are or are not.
   """
   def are_connected(self,idxAt1,idxAt2):
      sPair = '%d,%d' % (idxAt1,idxAt2)
      if self.dictPairs.has_key(sPair):
         sBond = self.dictPairs[sPair]
         if sBond in ['-','=','#','$',':','~']:
            return 1
      return 0
   
   """------------------------------------------------------------------ 
       shortestPaths: find all shortest paths between atoms idxAt1 and
                      idxAt2. The number of shortest paths is lower or
                      equal to nRings+1: an acyclic molecule has exactly
                      one shortest path between any two atoms.

       return: lstPaths = [lst1,lst2,...], which is a list of shortest
                          paths between atoms idxAt1 and idxAt2
                          lst = sequence of indexes for atoms starting
                          with idxAt1 neighbor towards idxAt2;
                          lst is empty if atoms idxAt1 and idxAt2 are
                          adjacent (or not at all connected via a path)
               None, if atom indexes out of range         
   """
   def shortest_paths(self,idxAt1,idxAt2):
      lstPaths = []
      nDist = self.dist(idxAt1,idxAt2)
      if nDist < 2:
         return lstPaths # return empty list if adjacent
      else:
         """
             start with all adjacent atoms, for which the distance to
             atom idxAt2 is less that the idxAt1-idxAt2 distance
         """
         if self.dictNbors.has_key(idxAt1):
            lstNbors =  self.dictNbors[idxAt1]
            for iNbor in lstNbors:
               if self.dist(iNbor,idxAt2) < nDist: 
                  lstGrowingPath =  [iNbor]
                  lstPaths.append(lstGrowingPath)
         else:
            return None
      """
         For distance nDist between atoms idxAt1 and idxAt2 the number
         of atoms in the shortest-path-sequence is nDist-1.
      """
      for k in range(nDist-1):
         lstPathsTemp = []
         for lstGrowingPath in lstPaths:
            lenPath = len(lstGrowingPath) 
            iEnd = lstGrowingPath[lenPath-1]
            if self.dictNbors.has_key(iEnd):
               lstEndNbors = self.dictNbors[iEnd]
               for iEndNbor in lstEndNbors:
                  if iEndNbor == iEnd:     # skip backward path
                     continue
                  elif iEndNbor == idxAt2: # obtain shortest path
                     lstPathsTemp.append(lstGrowingPath)
                  elif self.dist(iEndNbor,idxAt2) > nDist-(k+2): 
                     continue              # reject any longer path
                  else:                   
                     if k < nDist-2:       # still growing
                        lst = lstGrowingPath[:]
                        lst.append(iEndNbor)
                        lstPathsTemp.append(lst)
            else:
               # FATAL ERROR: bail out
               sMsg  = 'csm_subnotation.shortestPaths: atom with index ' 
               sMsg += '%d not in dictNbors while searching %d-%d paths'\
                       % (iEnd,idxAt1,idxAt2)
               self.lstErrors.append(sMsg)             
               return None
         lstPaths = lstPathsTemp
      return lstPaths

   #===================================================================#
   # ACCESS members (instead of direct access)                         #
   #===================================================================#
   def get_sAnnSmi(self):       return self.sAnnSmi
   def get_lstTokens(self):     return self.lstTokens 
   def get_lstTokTyp(self):     return self.lstTokTyp 
   def get_lstDepth(self):      return self.lstDepth  
   def get_maxDepth(self):      return self.maxDepth  
   def get_nNodes(self):        return self.nNodes
   def get_nAtoms(self):        return self.nAtoms
   def get_nHterm(self):        return self.nHterm 
   def get_nAhold(self):        return self.nAhold     
   def get_nRings(self):        return self.nRings
   def get_nLocal(self):        return self.nLocal    
   def get_nDeloc(self):        return self.nDeloc 
   def get_lstAtSymb(self):     return self.lstAtSymb  
   def get_lstAtNumb(self):     return self.lstAtNumb   
   def get_lstAtLabel(self):    return self.lstAtLabel 
   def get_lstAtCharge(self):   return self.lstAtCharge  
   def get_lstAtDepth(self):    return self.lstAtDepth     
   def get_lstAromat(self):     return self.lstAromat
   def get_lstAaaEntr(self):    return self.lstAaaEntr
   def get_lstHcount(self):     return self.lstHcount
   def get_lstNbors(self):      return self.lstNbors
   def get_lstNvbors(self):     return self.lstNvbors
   def get_lstRingMap(self):    return self.lstRingMap  
   def get_lstCaaEntr(self):    return self.lstCaaEntr
   def get_lstAtf0(self):       return self.lstAtf0      
   def get_dictNbors(self):     return self.dictNbors 
   def get_dictBonds(self):     return self.dictBonds 
   def get_dictPairs(self):     return self.dictPairs 
   def get_lstRingLength(self): return self.lstRingLength 
   def get_lstRingAromat(self): return self.lstRingAromat
   def get_lstRingCharge(self): return self.lstRingCharge
   def get_lstRings(self):      return self.lstRings
   def get_lstIdxMat(self):     return self.lstIdxMat
   def get_lstAdjMat(self):     return self.lstAdjMat
   def get_lstDistMat(self):    return self.lstDistMat 
   def get_lstErrors(self):     return self.lstErrors

   #===================================================================#
   # ACCESS molecular-graph data by atom list methods                  #
   #===================================================================#
   def atf0(self):          return self.lstAtf0
   def atoms_aromat(self):  return self.lstAromat
   def atoms_charge(self):  return self.lstAtCharge
   def atoms_label(self):   return self.lstAtLabel    
   def atoms_hcount(self):  return self.lstHcount  
   def atoms_nbors(self):   return self.lstNbors
   def atoms_nvbors(self):  return self.lstNvbors   
   def atoms_numb(self):    return self.lstAtNumb   
   def atoms_symb(self):    return self.lstAtSymb
   def atoms_nvbors(self):  return self.lstNvbors   


   #===================================================================#
   # ACCESS molecular-graph data by component methods                  #
   #===================================================================#
   """------------------------------------------------------------------   
      aaa_entries:

      EXAMPLE: self.AnnSmi = 'C1COCC[N+]{+Rn=1,2}{+R}';
               (encoding 4-methyl- and 4-ethyl-4-alkylmorpholinium)
               iuAtom = 6
               returns [ ['+R',{'n':'1,2'} ], ['+R',{}] ] f

     return: [ [sAM1,dictAnn1],  ... ] 
   """
   def aaa_entries(self,iuAtom):
      if iuAtom > 0 and iuAtom <= self.nNodes:
         lstObj =  self.lstAaaEntr[iuAtom-1]
         lstEntries = []
         for oCurly in lstObj:
            lstEntries.append(oCurly.entry())
         return lstEntries
      else:
         return None

   """------------------------------------------------------------------   
      aaa_indices:
      return: list with indices of those atoms to which annotations are
              anchored
              (indices between (including) 0 and self.nAtoms-1)              
   """   
   def aaa_indices(self):
      lstIndices = []
      iAtom = 0
      for lst in self.lstAaaEntr:
         if len(lst) > 0:
            lstIndices.append(iAtom)
         iAtom += 1
      return lstIndices

   """------------------------------------------------------------------   
      aaa_indices_iu:
      return: list with indices of those atoms to which annotations are
              anchored
              (indices between (including) 1 and self.nAtoms)
   """   
   def aaa_indices_iu(self):
      lstIndices = []
      iAtom = 0
      for lst in self.lstAaaEntr:
         if len(lst) > 0:
            lstIndices.append(iAtom+1)
         iAtom += 1
      return lstIndices

   def caa_entries(self):
      lstEntries = []
      for oCurly in self.lstCaaEntr:
         lstEntries.append(oCurly.entry())
      return lstEntries
   

   """------------------------------------------------------------------   
      rings: get list of list with indices of ring members
                (indices between (including) 0 and self.nAtoms-1)
      return:   self.lstRings
   """   
   def rings(self):
      return self.lstRings


   """------------------------------------------------------------------   
      rings_iu: get list of list with indices of ring members
                (indices between (including) 1 and self.nAtoms)
      return:   self.lstRings_iu
   """   
   def rings_iu(self):
      lstRings_iu = []
      for lstMemb in self.lstRings:
         lstMemb_iu = []
         for iAtom in lstMemb:
            lstMemb_iu.append(iAtom+1)
         lstRings_iu.append(lstMemb_iu)
      return lstRings_iu

   """------------------------------------------------------------------   
      atpair_bond: get bond symbol for bond between atons iu and ju
                   (iu and ju between (including) 1 and self.nAtoms)
      return: sBond
   """   
   def atpair_bond(self,iu,ju):
      sPair = '%d,%d' % (iu-1,ju-1)
      if self.dictPairs.has_key(sPair):
         return self.dictPairs[sPair]
      else:
         return None
   
   def deloc_charge(self): return self.nDeloc 
   def dmat(self):         return self.lstDistMat
   def numof_atoms(self):  return self.nAtoms   
   def numof_nodes(self):  return self.nNodes
   def numof_rings(self):  return self.nRings   
   def msgs_err(self):     return self.lstErrors

   #===================================================================#
   # ACCESS MOLECULAR FORMULA data derived for current notation        #
   #===================================================================#
   def mf_linear_notation(self):    return self.oMf.linear_notation()
   def mf_hill_format(self):        return self.oMf.hill_format()   
   def mf_list_of_pairs(self):      return self.oMf.mf_list_of_pairs()
   def mf_dict_of_dict(self):       return self.oMf.dict_of_dict()
   def mf_charge_notation(self):    return self.oMf.charge_notation()
   def mf_charge_number(self):      return self.oMf.charge_number()
   def mf_msgs_err(self):           return self.oMf.msgs_err()
   
