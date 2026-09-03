# Makefile for the PoolParty resubmission and bioRxiv preprint.
#
LATEXMK   := latexmk -pdf -interaction=nonstopmode -halt-on-error
LATEXDIFF := latexdiff

RESUB_SUPPORT := references.bib sn-jnl.cls sn-vancouver-num.bst
BIORXIV_SUPPORT := ../legacy/26.04.06_manuscript/kinneylab_preprintstyle.cls \
	../legacy/26.04.06_manuscript/kinneylab_bibstyle.bst
ASSETS := $(wildcard assets/*)

PDFS := si_resub.pdf main_resub.pdf si_biorxiv.pdf main_biorxiv.pdf response.pdf
CLEANPDFS := si_resub-clean.pdf main_resub-clean.pdf
NOCOLOUR := -pretex='\def\NOQUOTEDCOLOR{}' -usepretex

# Log messages meaning a PDF is not submission-ready.
BADLOG := 'There were undefined (references|citations)|Reference .* undefined|Citation .* undefined|multiply defined|Float too large|^! |Emergency stop|Fatal error occurred'

.PHONY: all si resub main biorxiv response tracked clean-pdf verify clean distclean help
.DEFAULT_GOAL := all

all: $(PDFS) verify        ## Build all resubmission and bioRxiv PDFs

si: si_resub.pdf           ## Build the resubmission Supplemental Information
resub: main_resub.pdf      ## Build the resubmission manuscript (and SI)
main: main_resub.pdf       ## Alias for `make resub`
biorxiv: main_biorxiv.pdf ## Build the bioRxiv manuscript (and SI)
response: response.pdf     ## Build the response letter (and resubmission documents)
tracked: main_tracked.pdf  ## Build the tracked-changes resubmission manuscript

si_resub.pdf: si_resub.tex si_body.tex $(RESUB_SUPPORT) $(ASSETS)
	$(LATEXMK) si_resub.tex

main_resub.pdf: main_resub.tex main_body.tex $(RESUB_SUPPORT) $(ASSETS)
	$(LATEXMK) main_resub.tex

response.pdf: response.tex $(RESUB_SUPPORT)
	$(LATEXMK) response.tex

# bioRxiv chain. .latexmkrc makes the legacy class and bibliography style
# available to `latexmk`, including editor-triggered builds.
si_biorxiv.pdf: si_biorxiv.tex si_body.tex $(RESUB_SUPPORT) $(ASSETS)
	$(LATEXMK) si_biorxiv.tex

main_biorxiv.pdf: main_biorxiv.tex main_body.tex references.bib $(BIORXIV_SUPPORT) $(ASSETS)
	$(LATEXMK) main_biorxiv.tex

# Unmarked resubmission copies for submission.
clean-pdf: $(CLEANPDFS)   ## Build unmarked resubmission manuscript and SI PDFs

si_resub-clean.pdf: si_resub.tex si_body.tex $(RESUB_SUPPORT) $(ASSETS)
	$(LATEXMK) $(NOCOLOUR) -jobname=si_resub-clean si_resub.tex

main_resub-clean.pdf: main_resub.tex main_body.tex $(RESUB_SUPPORT) $(ASSETS)
	$(LATEXMK) $(NOCOLOUR) -jobname=main_resub-clean main_resub.tex

# Tracked-changes manuscript, diffed against the originally submitted version.
main_tracked.tex: baseline_submission.tex main_resub.tex main_body.tex
	$(LATEXDIFF) baseline_submission.tex main_resub.tex > $@.tmp && mv $@.tmp $@

main_tracked.pdf: main_tracked.tex main_body.tex $(RESUB_SUPPORT) $(ASSETS)
	$(LATEXMK) main_tracked.tex

verify:                    ## Check logs for unresolved references/citations and fatal layout errors
	@fail=0; for f in si_resub main_resub si_biorxiv main_biorxiv response main_tracked si_resub-clean main_resub-clean; do if [ -f $$f.log ]; then if grep -Eq $(BADLOG) $$f.log; then echo "FAIL  $$f.log"; grep -En $(BADLOG) $$f.log; fail=1; else echo "  ok  $$f.log"; fi; fi; done; exit $$fail

clean:                     ## Remove LaTeX auxiliary files, keeping PDFs
	$(LATEXMK) -bibtex -c si_resub.tex main_resub.tex si_biorxiv.tex main_biorxiv.tex response.tex
	@rm -f main_tracked.tex main_tracked.tex.tmp main_tracked.aux main_tracked.log main_tracked.out main_tracked.fls main_tracked.fdb_latexmk main_tracked.bbl main_tracked.blg main_tracked.toc
	@rm -f si_resub-clean.aux si_resub-clean.log si_resub-clean.out si_resub-clean.fls si_resub-clean.fdb_latexmk si_resub-clean.bbl si_resub-clean.blg si_resub-clean.toc
	@rm -f main_resub-clean.aux main_resub-clean.log main_resub-clean.out main_resub-clean.fls main_resub-clean.fdb_latexmk main_resub-clean.bbl main_resub-clean.blg

distclean: clean           ## Also delete generated PDFs
	$(LATEXMK) -bibtex -C si_resub.tex main_resub.tex si_biorxiv.tex main_biorxiv.tex response.tex
	@rm -f main_tracked.pdf $(CLEANPDFS)

help:                      ## Show available build targets
	@grep -E '^[a-zA-Z_.-]+:.*## ' $(MAKEFILE_LIST) | sed 's/:[^#]*## /|/' | awk -F'|' '{printf "  make %-10s %s\n", $$1, $$2}'
