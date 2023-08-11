simulated_chrom_21.ts.tsz:
	wget https://zenodo.org/record/7702392/files/simulated_chrom_21.ts.tsz

# This is *really* big
simulated_chrom_1.ts.tsz:
	wget https://zenodo.org/record/7702392/files/simulated_chrom_1.ts.tsz

%.ts: %.ts.tsz
	tsunzip -k $<

%.bcf: %.ts
	python3 -m tskit vcf $< | bcftools view -O b > $@

%.sgz: %.bcf
	python3 ../src/convert.py vcf-to-sgkit-zarr $< $@
