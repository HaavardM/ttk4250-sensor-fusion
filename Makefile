


run:
	matlab -nosplash -nodisplay -r "run('matlab/main.m');exit;" | tail -n +11

build-pdf:
	mkdir -p latex/build
	cd latex  && pdflatex --output-directory ./build ov.tex
	cp latex/build/ov.pdf results.pdf
clean:
	rm -rf latex/build

.PHONY: clean
