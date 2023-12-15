BCFTOOLS=../software/bcftools
BCFTOOLS_PLUGINS=../software/bcftools-1.18/plugins
SAVVY=../software/savvy/bin/sav
GENOZIP=../software/genozip

simulated_chrom_21.ts.tsz:
	wget https://zenodo.org/record/7702392/files/simulated_chrom_21.ts.tsz

# This is *really* big
simulated_chrom_1.ts.tsz:
	wget https://zenodo.org/record/7702392/files/simulated_chrom_1.ts.tsz

%.ts: %.ts.tsz
	tsunzip -k $<

.PRECIOUS: data/%.bcf
%.bcf: %.ts
	python3 -m tskit vcf $< | ${BCFTOOLS} view -O b > $@
	${BCFTOOLS} index $@

# Compute some minimal INFO summaries to support the af-dist computation
%.tags.bcf: %.bcf
	BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS} ${BCFTOOLS} \
			 +fill-tags $< -Ob -o $@ -- -t AN,AC,AF
	${BCFTOOLS} index $@

%.sgz: %.bcf
	python3 ../src/convert.py vcf-to-sgkit-zarr $< $@

%.sav: %.bcf
	${SAVVY} import --phasing full $< > $@

%.genozip: %.bcf
	# There's some weird bug in genocat that doesn't allow us 
	# to run --samples on output if it's derived directly from
	# bcf, so we have to filter it through bcftools first.
	${BCFTOOLS} view $< | ${GENOZIP} --no-test -o $@

# NOTES here for now
# it's quite hard to find something that bcftools does that just
# outputs one thing. Looking through the bcftools plugins, it
# seems that af-dist might be something that we can compare to
# without doing something totally apples-to-oranges
#
# This might actually be a good opportunity to contrast the 
# development model of a C API with numba etc.
# BCFTOOLS_PLUGINS=../software/bcftools-1.18/plugins ../software/bcftools +af-dist
#
%.bcf.sstats: %.bcf
	# /usr/bin/time
	${BCFTOOLS} stats -s - $< > $@
