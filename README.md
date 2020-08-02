# snpview

*INSTALL*
sudo apt-get install libncurses5-dev
cd samtools-0.1.19
make
cd ..
make

*USAGE*
```
name:	snpview 0.2
func:	a text based multiple BAM alignment viewer
author:	Yi Wang @ Fudan University
usage:	snpview [options] <chr pos> [bam1 bam2...]
	-l <bams>	add BAMs listed in a file
	-r <ref.fa>	load reference
	-w <width>	half window size
	-R <val>	set max base quality in red (3)
	-Y <val>	set max base quality in yellow (6)
	-G <val>	set max base quality in green (18)
	-c		use vt100 color
	-C		collapse insertions
	-v		show variant reads only
	-n		don't reset the console
	-m		filter read by mask
```
