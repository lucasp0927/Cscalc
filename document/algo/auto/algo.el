(TeX-add-style-hook "algo"
 (lambda ()
    (LaTeX-add-labels
     "eq:timeevo"
     "eq:motion"
     "fig:pulses"
     "fig:pulse"
     "eq:perturb"
     "eq:f")
    (TeX-run-style-hooks
     "lmodern"
     "listings"
     "geometry"
     "url"
     "subfig"
     "float"
     "graphicx"
     "graphics"
     "xeCJK"
     "bm"
     "amssymb"
     "amsmath"
     "fontspec"
     "latex2e"
     "art11"
     "article"
     "11pt"
     "a4paper")))

