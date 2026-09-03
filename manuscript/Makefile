# Makefile for the PoolParty BMC Bioinformatics submission.
#
# Build order is significant. main.tex imports si.tex's labels via xr-hyper,
# and response.tex imports labels from both main.tex and si.tex, so each
# document reads the .aux file of the ones before it:
#
#     si.tex  -->  main.tex  -->  response.tex
#
# The rules below encode that chain, so editing si.tex correctly forces
# main.pdf and response.pdf to rebuild with the new Figure/Table S* numbers.
# .NOTPARALLEL keeps the ordering intact even under "make -j".
#
# Run "make help" for the available targets.

LATEXMK   := latexmk -pdf -interaction=nonstopmode -halt-on-error
LATEXDIFF := latexdiff

# Shared inputs: a change to any of these rebuilds every document.
SUPPORT := references.bib sn-jnl.cls sn-vancouver-num.bst
ASSETS  := $(wildcard assets/*)

PDFS := si.pdf main.pdf response.pdf

# Unmarked copies for submission: same sources, quoted-text colouring off.
CLEANPDFS := si-clean.pdf main-clean.pdf
NOCOLOUR  := -pretex='\def\NOQUOTEDCOLOR{}' -usepretex

# Log messages meaning the PDF built but is not submission-ready.
BADLOG := 'There were undefined (references|citations)|Reference .* undefined|Citation .* undefined|multiply defined|Float too large'

.PHONY: all si main response tracked clean-pdf verify clean distclean help
.NOTPARALLEL:
.DEFAULT_GOAL := all

all: $(PDFS) verify        ## Build all three PDFs, then check their logs

si: si.pdf                 ## Build the Supplemental Information only
main: main.pdf             ## Build the manuscript (builds si.pdf first)
response: response.pdf     ## Build the response letter (builds si.pdf, main.pdf first)
tracked: main_tracked.pdf  ## Build the latexdiff tracked-changes manuscript

# Supplemental Information: no external document dependencies.
si.pdf: si.tex $(SUPPORT) $(ASSETS)
	$(LATEXMK) si.tex

# Manuscript: reads si.aux for the Figure S* and Table S* numbers.
main.pdf: main.tex si.pdf $(SUPPORT) $(ASSETS)
	$(LATEXMK) main.tex

# Response letter: reads main.aux and si.aux so quoted cross-references resolve.
response.pdf: response.tex main.pdf si.pdf $(SUPPORT)
	$(LATEXMK) response.tex

clean-pdf: $(CLEANPDFS)   ## Build unmarked si-clean.pdf/main-clean.pdf for submission

# main.tex still reads si.aux for the Figure/Table S* numbers, so plain
# si.pdf stays a prerequisite; those numbers do not depend on colouring.
main-clean.pdf: main.tex si.pdf $(SUPPORT) $(ASSETS)
	$(LATEXMK) $(NOCOLOUR) -jobname=main-clean main.tex

si-clean.pdf: si.tex $(SUPPORT) $(ASSETS)
	$(LATEXMK) $(NOCOLOUR) -jobname=si-clean si.tex

# Tracked-changes manuscript, diffed against the originally submitted version.
main_tracked.tex: baseline_submission.tex main.tex
	$(LATEXDIFF) baseline_submission.tex main.tex > $@.tmp && mv $@.tmp $@

main_tracked.pdf: main_tracked.tex si.pdf $(SUPPORT) $(ASSETS)
	$(LATEXMK) main_tracked.tex

verify:                    ## Check build logs for undefined refs/citations, bad floats
	@fail=0; for f in si main response main_tracked si-clean main-clean; do if [ -f $$f.log ]; then if grep -Eq $(BADLOG) $$f.log; then echo "FAIL  $$f.log"; grep -En $(BADLOG) $$f.log; fail=1; else echo "  ok  $$f.log"; fi; fi; done; exit $$fail

clean:                     ## Remove LaTeX auxiliary files, keeping the PDFs
	$(LATEXMK) -c si.tex main.tex response.tex
	@rm -f main_tracked.tex main_tracked.tex.tmp main_tracked.aux main_tracked.log main_tracked.out main_tracked.fls main_tracked.fdb_latexmk main_tracked.bbl main_tracked.blg main_tracked.toc
	@rm -f si-clean.aux si-clean.log si-clean.out si-clean.fls si-clean.fdb_latexmk si-clean.bbl si-clean.blg si-clean.toc main-clean.aux main-clean.log main-clean.out main-clean.fls main-clean.fdb_latexmk main-clean.bbl main-clean.blg

distclean: clean           ## Also delete the PDFs (tracked in git; restore with git checkout)
	$(LATEXMK) -C si.tex main.tex response.tex
	@rm -f main_tracked.pdf $(CLEANPDFS)

help:                      ## Show this help
	@grep -E '^[a-zA-Z_.-]+:.*## ' $(MAKEFILE_LIST) | sed 's/:[^#]*## /|/' | awk -F'|' '{printf "  make %-10s %s\n", $$1, $$2}'
