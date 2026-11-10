# GitHub keeps TeX sources and exact proof data, not generated PDFs.
# Run from the repository root: latexmk papers/<manuscript>.tex
# Also works with 'latexmk -pdf' (we still invoke LuaLaTeX).
$pdf_mode = 4;
$lualatex = 'lualatex -interaction=nonstopmode -halt-on-error -file-line-error -synctex=0 %O %S';
$pdflatex = $lualatex;
$max_repeat = 6;
$out_dir = 'work/latexmk';
