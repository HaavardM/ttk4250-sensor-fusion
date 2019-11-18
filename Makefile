


task1:
	matlab -nosplash -nodisplay -r "try, run('task1/main.m'), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" | tail -n +11
	cp -r task1/plots/* latex/plots/a1/
task2:
	matlab -nosplash -nodisplay -r "try, run('task2/main.m'), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" | tail -n +11
build-pdf: task1 task2
	mkdir -p latex/build
	cd latex  && pdflatex --output-directory ./build ov.tex
	cd latex  && pdflatex --output-directory ./build ov.tex
	cp latex/build/ov.pdf results.pdf
clean:
	rm -rf latex/build

.PHONY: clean
.PHONY: task1
.PHONY: task2
