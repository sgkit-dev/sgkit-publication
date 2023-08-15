bcftools:
	wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
	tar -jxvf bcftools-1.18.tar.bz2
	cd bcftools-1.18 && ./configure && make -j 8
	cp bcftools-1.18/bcftools ./


clean:
	rm -fR bcftools*
