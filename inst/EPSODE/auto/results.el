(TeX-add-style-hook
 "results"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8" "latin1") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "listings"
    "color"
    "amsmath"
    "array"
    "fontenc"
    "natbib"
    "authblk"
    "inputenc"
    "amssymb"
    "booktabs"
    "blindtext"
    "scrextend"
    "geometry"
    "graphicx"
    "grffile"
    "longtable"
    "wrapfig"
    "rotating"
    "ulem"
    "textcomp"
    "capt-of"
    "hyperref"))
 :latex)

