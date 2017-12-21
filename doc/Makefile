TeX=xelatex
latexfile=doc

$(latexfile).pdf: *.tex
	while ($(TeX) $(latexfile); \
		grep -q "Rerun to get cross-references right." $(latexfile).log ) do true;\
	done

all: $(latexfile)
$(latexfile): *.tex
	$(TeX) $(latexfile);
	#bibtex $(latexfile);
	#$(TeX) $(latexfile);
	$(TeX) $(latexfile);

clean:
	rm -f *.log *.out *.blg *.aux *.toc $(latexfile).dvi $(latexfile).ps
