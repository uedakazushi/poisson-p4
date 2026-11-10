# Shared latexmk settings for CLI and VS Code LaTeX Workshop.
# The source files use \input{papers/...} relative to the repository root.
# Always run latexmk with the repository root as the working directory.
$pdf_mode = 4;
$lualatex = 'lualatex -interaction=nonstopmode -halt-on-error -file-line-error -synctex=1 %O %S';
$pdflatex = $lualatex;
$max_repeat = 6;
$out_dir = 'work/latexmk';
$do_cd = 0;
use File::Path qw(make_path);
make_path($out_dir) unless -d $out_dir;
