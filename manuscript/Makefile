# Makefile for the PoolParty BMC Bioinformatics submission.
#
# Build order is significant. main_bmc.tex imports si.tex's labels via
# xr-hyper, and response.tex imports labels from both the BMC manuscript and
# si.tex, so each
# document reads the .aux file of the ones before it:
#
#     si.tex  -->  main_bmc.tex  -->  response.tex
#                 main_biorxiv.tex
#
# The rules below encode that chain, so editing si.tex correctly forces
# main.pdf, main_biorxiv.pdf, and response.pdf to rebuild with the new
# Figure/Table S* numbers.
# .NOTPARALLEL keeps the ordering intact even under "make -j".
#
# Run "make help" for the available targets.

LATEXMK   := latexmk -pdf -interaction=nonstopmode -halt-on-error
LATEXDIFF := latexdiff

# Shared inputs: a change to any of these rebuilds every document.
SUPPORT := references.bib sn-jnl.cls sn-vancouver-num.bst
BIORXIV_SUPPORT := ../legacy/26.04.06_manuscript/kinneylab_preprintstyle.cls \
	../legacy/26.04.06_manuscript/kinneylab_bibstyle.bst
# TEXINPUTS/BSTINPUTS for the two files above now live in .latexmkrc, so an
# editor-triggered build finds them too, not just one launched from here.
ASSETS  := $(wildcard assets/*)

PDFS := si.pdf main.pdf main_biorxiv.pdf response.pdf

# Unmarked copies for submission: same sources, quoted-text colouring off.
CLEANPDFS := si-clean.pdf main-clean.pdf
NOCOLOUR  := -pretex='\def\NOQUOTEDCOLOR{}' -usepretex

# Log messages meaning the build hard-failed, or the PDF built but is not
# submission-ready. The hard-failure patterns matter: without them "verify"
# reported "ok" for a log whose run ended in "Fatal error occurred, no output
# PDF file produced!", green-lighting a document that never built.
BADLOG := 'There were undefined (references|citations)|Reference .* undefined|Citation .* undefined|multiply defined|Float too large|^! |Emergency stop|Fatal error occurred'

.PHONY: all si main biorxiv response tracked clean-pdf verify clean distclean help
.NOTPARALLEL:
.DEFAULT_GOAL := all

all: $(PDFS) verify        ## Build all four PDFs, then check their logs

si: si.pdf                 ## Build the Supplemental Information only
main: main.pdf             ## Build the BMC manuscript (builds si.pdf first)
biorxiv: main_biorxiv.pdf  ## Build the unmarked bioRxiv manuscript (builds si.pdf first)
response: response.pdf     ## Build the response letter (builds si.pdf, main.pdf first)
tracked: main_tracked.pdf  ## Build the latexdiff tracked-changes manuscript

# Supplemental Information: no external document dependencies.
si.pdf: si.tex $(SUPPORT) $(ASSETS)
	$(LATEXMK) si.tex

# BMC manuscript: preserve the established main.pdf/main.aux names so the
# response letter and existing submission workflow continue to work.
main.pdf: main_bmc.tex main_body.tex si.pdf $(SUPPORT) $(ASSETS)
	$(LATEXMK) -jobname=main main_bmc.tex

# bioRxiv manuscript: the wrapper uses the legacy preprint class directly.
main_biorxiv.pdf: main_biorxiv.tex main_body.tex si.pdf references.bib $(BIORXIV_SUPPORT) $(ASSETS)
	$(LATEXMK) main_biorxiv.tex

# Response letter: reads main.aux and si.aux so quoted cross-references resolve.
response.pdf: response.tex main.pdf si.pdf $(SUPPORT)
	$(LATEXMK) response.tex

clean-pdf: $(CLEANPDFS)   ## Build unmarked si-clean.pdf/main-clean.pdf for submission

# main_bmc.tex still reads si.aux for the Figure/Table S* numbers, so plain
# si.pdf stays a prerequisite; those numbers do not depend on colouring.
main-clean.pdf: main_bmc.tex main_body.tex si.pdf $(SUPPORT) $(ASSETS)
	$(LATEXMK) $(NOCOLOUR) -jobname=main-clean main_bmc.tex

si-clean.pdf: si.tex $(SUPPORT) $(ASSETS)
	$(LATEXMK) $(NOCOLOUR) -jobname=si-clean si.tex

# Tracked-changes manuscript, diffed against the originally submitted version.
main_tracked.tex: baseline_submission.tex main_bmc.tex
	$(LATEXDIFF) baseline_submission.tex main_bmc.tex > $@.tmp && mv $@.tmp $@

main_tracked.pdf: main_tracked.tex si.pdf $(SUPPORT) $(ASSETS)
	$(LATEXMK) main_tracked.tex

verify:                    ## Check build logs for undefined refs/citations, bad floats
	@fail=0; for f in si main main_biorxiv response main_tracked si-clean main-clean; do if [ -f $$f.log ]; then if grep -Eq $(BADLOG) $$f.log; then echo "FAIL  $$f.log"; grep -En $(BADLOG) $$f.log; fail=1; else echo "  ok  $$f.log"; fi; fi; done; exit $$fail

# latexmk protects .bbl from cleanup unless bibtex_use=2, which "-bibtex"
# sets: a .bbl can be a document's only bibliography source. Here it is
# always regenerable from the tracked references.bib, so it is removed.
clean:                     ## Remove LaTeX auxiliary files, keeping the PDFs
	$(LATEXMK) -bibtex -c si.tex
	$(LATEXMK) -bibtex -c -jobname=main main_bmc.tex
	$(LATEXMK) -bibtex -c main_biorxiv.tex response.tex
	@rm -f main_tracked.tex main_tracked.tex.tmp main_tracked.aux main_tracked.log main_tracked.out main_tracked.fls main_tracked.fdb_latexmk main_tracked.bbl main_tracked.blg main_tracked.toc
	@rm -f si-clean.aux si-clean.log si-clean.out si-clean.fls si-clean.fdb_latexmk si-clean.bbl si-clean.blg si-clean.toc main-clean.aux main-clean.log main-clean.out main-clean.fls main-clean.fdb_latexmk main-clean.bbl main-clean.blg

distclean: clean           ## Also delete the PDFs (tracked in git; restore with git checkout)
	$(LATEXMK) -bibtex -C si.tex
	$(LATEXMK) -bibtex -C -jobname=main main_bmc.tex
	$(LATEXMK) -bibtex -C main_biorxiv.tex response.tex
	@rm -f main_tracked.pdf $(CLEANPDFS)

help:                      ## Show this help
	@grep -E '^[a-zA-Z_.-]+:.*## ' $(MAKEFILE_LIST) | sed 's/:[^#]*## /|/' | awk -F'|' '{printf "  make %-10s %s\n", $$1, $$2}'
