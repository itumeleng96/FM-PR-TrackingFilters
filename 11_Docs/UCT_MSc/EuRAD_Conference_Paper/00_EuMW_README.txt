
CHANGES

20251014 EuMW Paper Template V6 Release
	- LaTeX only: use xurl package to manage long urls in paper text and references.
	- WORD and LaTeX: Section Abstract: Limit your abstract to 150 words in a single paragraph.
	- WORD and LaTeX: Section III.C.4: All markers next to an author name must be placed consecutively without comma or space among them (eg *^2).
	- WORD only: Section VI.C: WORD tables must not be used to house captions, equations or references.
	- LaTeX only: removed use of the biblatex package, and return to use of IEEEtran.bst for references.
	- WORD and LaTeX: minor changes in text of reference items due to use of IEEEtran.bst.

----------------------------------------
Information for EuMW Authors Using LaTeX
----------------------------------------

1)  The file
        EuMW_Paper_LaTeX_Template_V6.pdf
    contains all required format instructions for an EuMW paper. Authors must 
    read this PDF file so that they understand the required format of an EuMW paper.
    This PDF is created by compiling:
        EuMW_Paper_LaTeX_Template_V6.tex
    using pdfLaTeX, XeLaTeX, or LuaLaTeX.

2)  The file 
        EuMW_Paper_LaTeX_Template_V6.tex
    demonstrates how to format an EuMW paper correctly.  You should make a copy of
    this file and use it as the template for your paper, and type your own text into it.

    EuMW_Paper_LaTeX_Template_V6 has been configured to be compiled by pdfLaTeX
    or XeLaTeX or LuaLaTeX.  The file automatically detects which engine is being used
    and invokes the appropriate fonts in each case.  If you are using a different LaTeX
    engine, then you must invoke your own fonts to match the required format.

    pdfLaTeX and LuaLaTeX engines will use the Nimbus Roman fonts as main text fonts, while
    XeLaTeX will use the TeX Gyre Termes fonts.  These fonts are a close match to Times Roman
    used by the MS-Word version of the template.  Users of other LaTeX systems must ensure that 
    Nimbus Roman or TeX Gyre Termes or Times Roman fonts are used.

    pdfLaTeX, XeLaTeX, LuaLaTeX are available in the free TeXLive system described in
    (6) below.

3)  Authors should make a copy of
        EuMW_Paper_LaTeX_Template_V6.tex
    and then use the copy as a template by inserting their own text into it.

    Do not reuse your past .tex files as a template, even if your past papers conformed to 
    the required format.

4)  There is no requirement to balance both columns on any page including the last page, 
    i.e. both columns are not required to be vertically aligned at the bottom. 

    As shown in the file:
        EuMW_Paper_LaTeX_Template_V6.tex
    the following command must be used to prevent LaTeX from varying the vertical spacing
    to align the bottoms of both columns:
     \raggedbottom

5)  EuMW_Paper_LaTeX_Template_V6.tex uses the biblatex package to format the
    references section of the paper, by emulating the IEEE references style.

    The biblatex package is used instead of the common IEEEtran.bst file, because 
    biblatex offers better results when typesetting long URLs which are now common
    in reference lists.  You can see an example or a long URL in [1] of the template paper PDF.

6)  The template paper file 
        EuMW_Paper_LaTeX_Template_V6.tex 
    has been tested with the following 3 LaTeX implementations in TeXLive (6.1, 6.2, 6.3 below)
    TeXLive is freely available from: http://www.tug.org/texlive/

6.1) pdfLaTeX, installed as part of TeXLive

    The Computer Modern fonts as default main fonts for pdfLaTeX must not be used.

    The Nimbus Roman fonts are based on Times Roman and are installed as part of 
    TeXLive and must be used as main fonts instead of the Computer Modern fonts.

    In order to use the Nimbus Roman fonts installed with TeXLive, the template
    uses the commands:
        \usepackage{t1enc}
        \usepackage{times}

    To start writing your paper, do:
        copy EuMW_Paper_LaTeX_Template_V6.tex mypaper.tex

    To compile mypaper.tex to mypaper.pdf using pdfLaTeX, type the following commands:
        pdflatex mypaper
        bibtex mypaper
        pdflatex mypaper
        pdflatex mypaper

6.2) XeLaTeX, installed as part of TeXLive

    The Latin Modern Roman fonts as default main fonts for XeLaTeX must not be used.

    The TeX Gyre Termes fonts are based on Times Roman and are installed as part of
    TeXLive and must be used as main fonts instead of the Latin Modern Roman fonts.

    In order to use the TeX Gyre Termes fonts installed with TeXLive, the template
    uses the commands:
        \usepackage{fontspec}
        \setmainfont{TeX Gyre Termes}
   
    To start writing your paper, do:
        copy EuMW_Paper_LaTeX_Template_V6.tex mypaper.tex

    To compile mypaper.tex to mypaper.pdf using XeLaTeX, type the following commands:
        xelatex mypaper
        bibtex mypaper
        xelatex mypaper
        xelatex mypaper

6.3) LuaLaTeX, installed as part of TeXLive

    LuaLaTeX is able to use either the Nimbus Roman or TeX Gyre Termes fonts.  By default
    the TeX Gyre Termes is used in the template paper.

    To start writing your paper, do:
        copy EuMW_Paper_LaTeX_Template_V6.tex mypaper.tex

    To compile mypaper.tex to mypaper.pdf using LuaLaTeX, type the following commands:
        lualatex mypaper
        bibtex mypaper
        lualatex mypaper
        lualatex mypaper

----------------------------------------------------------

MANIFEST of this ZIP file.

1) EuMW_Paper_LaTeX_Template_V6.pdf
        This PDF contains format instructions for EuMW papers. 
        This PDF was created from the source file called
            EuMW_Paper_LaTeX_Template_V6.tex.

2) EuMW_Paper_LaTeX_Template_V6.tex
        Authors should write their papers by editing this
        source file and adding their own text.

3) The following accompanying files are included in the ZIP and are 
   required to build the sample paper EuMW_Paper_LaTeX_Template_V6.pdf.

        IEEEtran.cls
            IEEE paper style file 1.8b from CTAN
            You must keep this file in the same folder as your paper .tex file.

        EuMW_modify_IEEEtran_18b_CTAN_V6.tex
            \input this file after \documentclass{IEEEtran} to change various
            settings of IEEEtran into EuMW format.
            You must keep this file in the same folder as your paper .tex file.

        IEEEabrv.bib
            IEEE publication title abbreviations.
            You can discard this file if not needed.

        IEEEexample.bib
            IEEE example bibliography file
            You can insert your own bibliographic entries into this file or
            discard this file if you have your own .bib file.

        figures\*
            figure .eps files for the sample paper
            You can discard this folder.

        00_EuMW_README.txt
            this file

