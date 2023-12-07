all: bcftools savvy-afdist/savvy-afdist

bcftools:
	wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
	tar -jxvf bcftools-1.18.tar.bz2
	cd bcftools-1.18 && ./configure && make -j 8
	cp bcftools-1.18/bcftools ./

savvy:
	python3 -m venv env
	bash env/bin/activate
	python -m pip install cget
	cget install --prefix savvy statgen/savvy

savvy-afdist/savvy-afdist: savvy
	cmake -DCMAKE_TOOLCHAIN_FILE=../savvy/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release -B savvy-afdist savvy-afdist
	make -C savvy-afdist


clean:
	rm -fR bcftools* savvy*


    