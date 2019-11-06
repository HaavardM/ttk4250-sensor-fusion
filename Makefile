


run:
	matlab -nosplash -nodisplay -r "try, run('task1/main.m'), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" | tail -n +11
task2:
	matlab -nosplash -nodisplay -r "try, run('task2/main.m'), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" | tail -n +11
build-pdf:
	mkdir -p latex/build
	cd latex  && pdflatex --output-directory ./build ov.tex
	cp latex/build/ov.pdf results.pdf
clean:
	rm -rf latex/build

.PHONY: clean
