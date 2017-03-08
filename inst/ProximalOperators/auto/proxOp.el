(TeX-add-style-hook
 "proxOp"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("tcolorbox" "most") ("inputenc" "utf8") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "1_MathNotations"
    "article"
    "art10"
    "listings"
    "color"
    "amsmath"
    "array"
    "fontenc"
    "natbib"
    "ifthen"
    "xifthen"
    "xargs"
    "amssymb"
    "latexsym"
    "amsfonts"
    "dsfont"
    "empheq"
    "tcolorbox"
    "geometry"
    "inputenc"
    "graphicx"
    "grffile"
    "longtable"
    "wrapfig"
    "rotating"
    "ulem"
    "textcomp"
    "capt-of"
    "hyperref")
   (LaTeX-add-labels
    "sec:org3cb9764"
    "sec:org615ea4a"
    "sec:org7754336"
    "sec:orgc2478e1"
    "sec:orgfbb51e0"
    "sec:orgc1a56ca"
    "sec:org1b9b5bf"
    "sec:org6141ee3"
    "sec:orge51319d"
    "sec:orga180bfc"
    "sec:org65e550e"
    "sec:org02a150e"
    "sec:orgb73b055"
    "sec:orgbeacd7a"
    "sec:org539702e"
    "sec:org78ca640"
    "sec:org2b277c9")
   (LaTeX-add-bibliographies
    "biblio")
   (LaTeX-add-xcolor-definecolors
    "blue"
    "red"
    "0,.5,0"
    "gray"
    "white"))
 :latex)

