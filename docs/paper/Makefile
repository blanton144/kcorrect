NAME= paper

all: ${NAME}.ps

${NAME}.ps: ${NAME}.dvi 
	dvips -N0 -o ${NAME}.ps ${NAME}.dvi

${NAME}.dvi: ${NAME}.tex figures.tex bib_${NAME}.tex 
	latex ${NAME}
	latex ${NAME}
	latex ${NAME}

clean:
	rm ${NAME}.ps *.log *.dvi

fullclean:
	rm ${NAME}.ps *.log *.dvi *.aux

aj: dummy
	${HOME}/environment/bin/mkaj

apj: dummy
	${HOME}/environment/bin/mkapj

common: dummy
	${HOME}/environment/bin/mkcommon

dummy:
