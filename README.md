****************************************************************************
Python modules CurlySMILES in Python |
Version 1.0.1 | Axeleratio's CurlySMILES Project: http://www.axeleratio.com/csm/proj/main.htm |
CurlySMILES publication: http://www.jcheminf.com/content/3/1/1
****************************************************************************

Date: November 9, 2010 (updated March 29, 2013)

0. License: GNU General Public License (see: COPYING)

I. Description

    This package of Python modules includes code to parse, interpret
    and integrate CurlySMILES notations into Python-based software
    including CGIs. 
    CurlySMILES is chemical annotation and query language,
    documented at www.axeleratio.com/csm/proj/language.htm
    and explained, by example, at 
    www.axeleratio.com/csm/proj/examples.htm
    
    For a first-glance introduction, run:  python code_snippet.py
                                           ----------------------
    ( also see: http://www.axeleratio.com/csm/proj/py/doc/codesnippet.pdf )

    Run python testing.py to test a diverse set of 85 notations.
        -----------------
    
II. List of module files

    o   csm_aliases.py  : contains class AliasNotations to verify aliases
                          and retrieve associated notations.
    o   csm_annsmi.py   : contains class AnnotatedSmiles to manage and
                          parse a SMILES notation and its curly-braces-
                          enclosed annotations.
    o   csm_atoms.py    : contains class AtomicData to identify atomic
                          symbols and to provide atomic data.
    o   csm_curlyann.py : contains class CurlyAnnotation to manage and
                          parse a curly-braces-enclosed annotation. 
    o   csm_dataface.py : contains class DataFace to access any data required 
                          for the interpretation of a CurlySMILES notation
                          (to be used as singleton).
    o   csm_molform.py  : contains class MolecularFormula to manage, parse
                          and build molecular formulae in linear notation 
                          format (including charge notation and isotopically 
                          labelled atomic symbols).    
    o   csm_notation.py : contains the core class Notation to manage and 
                          parse a CurlySMILES notation.
    o   csm_sfn.py      : contains class StoichFormNotation to manage and 
                          parse a stoichiometric formula notation (SFN).

    Further, three subdirectories with Python modules:

      o aliases : modules to lookup primary aliases, mainly for ions 
      o secalia : modules to lookup secondary aliases
      o keyval  : modules to verify values for selected keys in
                  annotation entries


III. Requirements

    1.  Python version 2.4.1 or higher

    2.  Operating system: platform independent (tested under Linux)


IV. Additional resources and notes

    For more information about the CurlySMILES language and its 
    applications go to:

    http://www.axeleratio.com/csm/proj/main.htm

    and 

    http://www.axeleratio.com/csm/proj/py/doc/CurlySMILESinPython.htm



    
