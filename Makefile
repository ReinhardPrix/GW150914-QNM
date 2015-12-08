SHELL = /bin/bash

.DELETE_ON_ERROR :

STATUS = $(STATUS_$(V))
STATUS_ = $(STATUS_0)
STATUS_0 = @echo "Making $@";

LATEX = ( latex -halt-on-error -file-line-error $* >/dev/null || { cat $*.log; exit 1; } )
PDFLATEX = ( pdflatex -halt-on-error -file-line-error $* >/dev/null || { cat $*.log; exit 1; } )
BIBTEX = ( bibtex $* >/dev/null || { cat $*.blg; exit 1; } )
DVIPS = dvips -Ppdf -q -t unknown $*.dvi >/dev/null
PS2PDF = ps2pdf -dAutoRotatePages=/None $*.ps $*.pdf >/dev/null
OCTAVE = octave -qfH

figsrcdir = ./figsources/
gnuplotcfg = $(wildcard $(figsrcdir)gnuplot.cfg)
extratex = git_tag.tex $(gnuplotcfg) $(generatedtables)


main = RingdownSearch.pdf
#
bibs = RingdownSearch.bib
#
texfigs = \
#
generatedtables =
#

.PHONY : main figs dist clean figclean

main : $(main)

figs :
	@$(MAKE) FIGS=1 main

tables :
	@$(MAKE) TABLES=1 main

clean :
	rm -f git_tag.tex $(main:.pdf=){.aux,.bbl,.blg,.log,.out,.pdf,*Notes.bib}

cleanfigs :
	rm -f $(texfigs) fig*.pdf && cd $(figsrcdir) &&  rm -f fig*-inc.eps && rm -f fig*.tex && rm -f fig*.pdf && rm -f fig*.dump

cleantables :
	rm -f $(generatedtables) && cd $(figsrcdir) && rm -f $(generatedtables) && rm -rf table*.dump

git_tag.tex ::
	$(STATUS)./git-tag.sh $@

cleandeep:
	cd $(figsrcdir) && rm -f *-inc.eps && rm -f *.tex && rm -f *.pdf && rm -f *.dump && rm -f $(generatedtables)

$(main) : %.pdf : %.tex $(extratex) $(bibs) $(texfigs) $(generatedtables)
	$(STATUS)$(PDFLATEX) && $(BIBTEX) && $(PDFLATEX) && $(PDFLATEX)

dist : $(main)
	$(STATUS)rm -rf dist/ && ../scripts/make_dist.pl $(main:.pdf=) dist/

## ----- generate LaTeX tables ----------
ifeq ($(TABLES),1)
$(generatedtables): $(addprefix $(figsrcdir), $(generatedtables))
	$(STATUS)cd $(figsrcdir) && cp $(generatedtables) ..
endif

## ----- generate actual figure plot files ----------
ifeq ($(FIGS),1)

$(texfigs) : %.pdf : $(figsrcdir)%.tex $(gnuplotcfg)
	$(STATUS)cd $(figsrcdir) && $(LATEX) && $(DVIPS) && $(PS2PDF) && rm -f $*.{aux,log,dvi,ps} && cp $*.pdf ..

$(texonlyfigs) : %.pdf : $(figsrcdir)%.tex $(gnuplotcfg)
	$(STATUS)cd $(figsrcdir) && $(LATEX) && $(DVIPS) && $(PS2PDF) && rm -f $*.{aux,log,dvi,ps} && cp $*.pdf ..

$(figsrcdir)fig_mis_%_CohMetric_LS.tex: $(figsrcdir)create_mismatch_plots_CohMetric_LS.m $(wildcard ./data/InjPars_CohMetric_LS*new_LogUnifEcc.dat)
	$(STATUS)cd $(figsrcdir) && $(OCTAVE) create_mismatch_plots_CohMetric_LS.m

$(figsrcdir)fig_mis_%_SemicohMetric_LS.tex: $(figsrcdir)create_mismatch_plots_SemicohMetric_LS.m $(wildcard ./data/InjPars_SemicohMetric_LS*new_LogUnifEcc.dat)
	$(STATUS)cd $(figsrcdir) && $(OCTAVE) create_mismatch_plots_SemicohMetric_LS.m

$(figsrcdir)fig_mis_%_CohMetric_SS.tex: $(figsrcdir)create_mismatch_plots_CohMetric_SS.m $(wildcard ./data/InjPars_CohMetric_SS*new_LogUnifEcc.dat)
	$(STATUS)cd $(figsrcdir) && $(OCTAVE) create_mismatch_plots_CohMetric_SS.m

$(figsrcdir)fig_mis_%_SemicohMetric_SS.tex: $(figsrcdir)create_mismatch_plots_SemicohMetric_SS.m $(wildcard ./data/InjPars_SemicohMetric_SS*new_LogUnifEcc.dat)
	$(STATUS)cd $(figsrcdir) && $(OCTAVE) create_mismatch_plots_SemicohMetric_SS.m

$(figsrcdir)fig_templates_per_dim_ScoX1.tex: $(figsrcdir)plot_templates_per_dim.m
	$(STATUS)cd $(figsrcdir) && $(OCTAVE) plot_templates_per_dim.m

$(figsrcdir)fig_max_Tsft.tex: $(figsrcdir)plot_max_Tsft.m
	$(STATUS)cd $(figsrcdir) && $(OCTAVE) plot_max_Tsft.m

$(figsrcdir)fig_ScoX1_requiredDepth.tex: $(figsrcdir)plot_ScoX1_requiredDepth.m
	$(STATUS)cd $(figsrcdir) && $(OCTAVE) plot_ScoX1_requiredDepth.m

$(figsrcdir)fig_uOccupancy_vs_ToP.tex: $(figsrcdir)plot_uOccupancy_vs_ToP.m
	$(STATUS)cd $(figsrcdir) && octapps_run plot_uOccupancy_vs_ToP.m

endif

## ---------- generate LaTeX table files ----------
