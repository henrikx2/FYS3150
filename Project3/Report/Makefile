BUILDDIR=build
FILENAME=report

pdf:
	#mkdir $(BUILDDIR) -p
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash \
	--to=latex \
	--output=$(FILENAME).pdf \
	--pdf-engine=xelatex \
	--toc \
	--variable urlcolor=cyan \
	--variable papersize=a4paper \
	--variable documentclass=article \
	--filter pandoc-citeproc \
	--bibliography=bibliography.bib \
	--csl=https://raw.githubusercontent.com/citation-style-language/styles/master/harvard-anglia-ruskin-university.csl \

