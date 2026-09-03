# Read by latexmk for every document in this directory, so the legacy bioRxiv
# class and bibliography style resolve no matter who invokes the build:
# "make", a bare "latexmk", or an editor integration such as VS Code's
# LaTeX Workshop. Previously only the Makefile exported these paths, so an
# editor-triggered build of main_biorxiv.tex died with
#   ! LaTeX Error: File `kinneylab_preprintstyle.cls' not found.
# and the resulting failure recorded in main_biorxiv.fdb_latexmk made every
# later "make" exit non-zero with "Nothing to do" until that file was deleted.
#
# ORDER IS LOAD-BEARING. ../legacy/26.04.06_manuscript holds a whole previous
# submission -- its own main.tex, main.aux, main.bbl and references.bib -- so
# putting it ahead of the working directory makes pdflatex read those stale
# files in place of this manuscript's. ensure_path prepends, so the '.' calls
# come second to land first in the path: working directory, then the legacy
# directory for the two files only it has, then the TeX tree (ensure_path
# appends the separator when the variable was unset, so kpathsea still
# searches its defaults).
ensure_path( 'TEXINPUTS', '../legacy/26.04.06_manuscript' );
ensure_path( 'BSTINPUTS', '../legacy/26.04.06_manuscript' );
ensure_path( 'TEXINPUTS', '.' );
ensure_path( 'BSTINPUTS', '.' );
