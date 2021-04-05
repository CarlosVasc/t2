# ---------------------------------------------------------
# type "make" command in Unix to create report.pdf 
# ---------------------------------------------------------

#list of LibreOffice Draw files
FILE=report
LODFIGS = $(patsubst %.odg,%.pdf,$(wildcard *.odg))

all: $(FILE).pdf
	evince $< &

$(FILE).pdf: $(wildcard *.tex) $(LODFIGS) python octave.log ngspice ngspice2 ngspice3.log ngspice4.log
	pdflatex  $(FILE).tex
	pdflatex  $(FILE).tex
	pdflatex  $(FILE).tex

python:
	python2 ../python/datag.py > testi.txt

octave.log:
	octave ../mat/rc.m > octave.log

ngspice: ngspice.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab1.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)

ngspice.log: ../sim/Ng1.net
	ngspice -b $< -o $@



ngspice2: ngspice2.log
	$(eval NGSFIGS2=$(shell grep _FIG2 $< | sed 's/_FIG2//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS2=$(addsuffix .ps, $(NGSFIGS2)))
	$(eval NGSTABS2=$(shell grep _TAB2 $< | sed 's/_TAB2//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach l, $(NGSTABS2), sed -n '/^$l_TAB2/,/^$l_END2/{p;/^$l_END2/q}' $< | grep -v $l_TAB2 | grep -v $l_END2 | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $l_tab2.tex;)
	$(foreach l, $(NGSFIGPS2), ps2pdf $l;)
	
ngspice2.log: ../sim/Ng2.cir
	ngspice -b $< -o $@
	

ngspice3.log: ../sim/Ng3.cir
	ngspice -b $< -o $@

ngspice4.log: ../sim/Ng4.cir
	ngspice -b $< -o $@

	
#convert libreoffice draw figures to pdf
%.pdf: %.odg
	soffice --headless --convert-to pdf $<




clean:
	@rm -rf *.aux *.bbl *.blg *.glg *.glo *.gls *.ilg *.ist *.lof *.txt
	@rm -rf *.log *.lot *.nlo *.nls *.out *.toc *~ *.*% ./*.pdf ./*.ps ./*.txt
	@rm -rf *_tab.tex octave-workspace *.eps *testi.txt *ngspice_1.tex *op_tab1.tex *op_tab2.tex *teorico1.tex *teorico2.tex
.PHONY: all clean ngspice octave


