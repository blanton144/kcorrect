NAME= paper

all: ${NAME}.ps

${NAME}.ps: ${NAME}.dvi ${HOME}/.redo
	dvips -N0 -o ${NAME}.ps ${NAME}.dvi

${NAME}.dvi: ${NAME}.tex bib_${NAME}.tex ${HOME}/.redo
	latex ${NAME}
	latex ${NAME}
	latex ${NAME}

bib_${NAME}.tex: ${NAME}.tex ${HOME}/environment/tex/biblioraw.tex \
									${HOME}/.redo
	auto_cite_nat ${NAME}.tex

cite: bib_${NAME}.tex

clean:
	rm ${NAME}.ps *.log *.dvi

fullclean:
	rm ${NAME}.ps *.log *.dvi *.aux

${HOME}/.redo:
	touch ${HOME}/.redo

touch:
	touch ${HOME}/.redo

again: touch ${NAME}.ps

aj: dummy
	${HOME}/environment/bin/mkaj

apj: dummy
	${HOME}/environment/bin/mkapj

common: dummy
	${HOME}/environment/bin/mkcommon

dummy:
