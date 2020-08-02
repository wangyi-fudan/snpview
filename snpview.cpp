/**
@file	snpview.cpp
@brief	view short reads over a given position in many BAMs
@author	Yi Wang
@date	06/29/2010
*/
#include	<algorithm>
#include	<dirent.h>
#include	<stdint.h>
#include	<unistd.h>
#include	<cstdlib>
#include	<cctype>
#include	<cstdio>
#include	<string>
#include	<vector>
#include	<sam.h>
#include	<zlib.h>
using	namespace	std;
typedef	unsigned	uint;
typedef	const	char	cchr;

class	SNPView{	
private:
	struct	Insertion{
		uint	pos,	len;
		vector<char>	str;
		vector<uint8_t>	qua;
	};
	uint	ref_count;
	bool	variant(const	bam1_t *b);	///<	whether a read carries a variant base at query position
	inline	void	art_char(char	Q,	char	R,	uint8_t	M);	///<	put Q to screen according to Reference and Mapping quality
public:
	bool	isc,	col,	var,	res;
	uint8_t	red,	yel,	gre;
	uint	min,	max,	pos,	wid,	msk;
	string	reg;
	vector<char>	ref;
	vector<string>	bam_file;
	
	static	int	fetch(const	bam1_t *b,	void *data);	///<	callback function for bam_fetch
	void	add_list(cchr	*F);
	void	load_ref(cchr	*F);
	void	document(void);
	void	options(int	ac,	char	**av);
	void	show_bam(cchr	*F);
};

bool	SNPView::variant(const	bam1_t *b){
	if(!ref.size()){	ref_count++;	return	true;	}
	uint	pseq=0,	pref=b->core.pos;
	uint8_t	*seq=bam1_seq(b);	
	uint32_t	*cp=bam1_cigar(b);
	uint16_t	cn=b->core.n_cigar;
	
	for(uint	c=0;	c<cn;	c++){
		uint	clen=(cp[c]>>BAM_CIGAR_SHIFT);
		switch(cp[c]&BAM_CIGAR_MASK){
		case	BAM_CMATCH:
			if(pos>=pref&&pos<pref+clen)	if(bam_nt16_rev_table[bam1_seqi(seq,	pseq+pos-pref)]!=toupper(ref[pos]))	return	true;
			pseq+=clen;	pref+=clen;	break;
		case	BAM_CINS:	case	BAM_CSOFT_CLIP:	
			if(pref==pos+1)	return	true;
			pseq+=clen;	break;
		case	BAM_CDEL:	case	BAM_CREF_SKIP:
			if(pos>=pref&&pos<pref+clen)	return	true;
			pref+=clen;	break;
		}
	}	
	ref_count++;
	return	false;
}

void	SNPView::art_char(char	Q,	char	R,	uint8_t	M){
	char	r=toupper(R);
	if(isc){
		if(M&128)	printf("\033[4m");
		if((M&127)>gre)	putchar(Q==r?'=':Q);
		else	if((M&127)>yel)	printf("\033[32m%c",	Q==r?'=':Q);
		else	if((M&127)>red)	printf("\033[33m%c",	Q==r?'=':Q);
		else	printf("\033[31m%c",	Q==r?'=':Q);
		printf("\033[0m");
	}
	else	putchar(Q==r?'=':Q);
}

int	SNPView::fetch(const	bam1_t *b,	void *data){
	SNPView	*me=(SNPView*)data;
	if(b->core.flag&me->msk)	return	0;
	if(me->var&&!me->variant(b))	return	0;
	
	bool	flagr=(me->ref.size()>0);
	cchr	*ref=flagr?&me->ref[0]:NULL;
	uint	s=me->min,	e=me->max,	p=me->pos;
	uint	pseq=0,	pref=b->core.pos;
	uint8_t	*seq=bam1_seq(b),	*qua=bam1_qual(b),	baseq=255;	
	uint32_t	*cp=bam1_cigar(b);
	uint16_t	cn=b->core.n_cigar;
	if(!cn)	return	0;
	
	vector<Insertion>	vins;
	Insertion	tins;
	
	for(uint	i=s;	i<pref;	i++)	putchar(' ');
	for(uint	c=0;	c<cn;	c++){
		uint	clen=(cp[c]>>BAM_CIGAR_SHIFT);
		switch(cp[c]&BAM_CIGAR_MASK){
		case	BAM_CMATCH:
			for(uint	i=0;	i<clen;	i++,	pseq++,	pref++)	if(pref>=s&&pref<=e){
				if(pref==p){	putchar('|');	baseq=qua[pseq];	}
				me->art_char(bam_nt16_rev_table[bam1_seqi(seq,	pseq)],	flagr?ref[pref]:'?',	qua[pseq]);
				if(pref==p)	putchar('|');
			}
			break;
		case	BAM_CINS:	case	BAM_CSOFT_CLIP:
			if(pref>=s+1&&pref<=e+1){
				tins.str.resize(clen);	tins.qua.resize(clen);	tins.pos=pref-1;	tins.len=clen;
				for(uint	i=0;	i<clen;	i++,	pseq++){
					tins.str[i]=bam_nt16_rev_table[bam1_seqi(seq,	pseq)];
					tins.qua[i]=qua[pseq];
				}
				vins.push_back(tins);
			}
			else	pseq+=clen;
			break;
		case	BAM_CDEL:	case	BAM_CREF_SKIP:
			for(uint	i=0;	i<clen;	i++,	pref++){
				if(pref==p)	printf("|*|");
				else	if(pref>=s&&pref<=e)	putchar('*');
			}
			break;
		}
	}
	for(uint	i=pref;	i<=e;	i++)	putchar(' ');
	
	if(baseq<255)	printf(" %2u",	baseq&127);
	else	printf(" NA");
	printf(" %3u%c\n",	b->core.qual,	b->core.flag&BAM_FREVERSE?'<':'>');
	
	if(!vins.size())	return	0;
	if(me->col){
		for(uint	j=s;	j<vins[0].pos;	j++)	printf("%s",	j==p?"   ":" ");
		for(uint	j=0;	j<vins.size();	j++){
			if(vins[j].pos==p)	putchar(' ');
			if(vins[j].len==1)	me->art_char(vins[j].str[0],	'?',	vins[j].qua[0]);
			else	if(vins[j].len<10)	putchar('0'+vins[j].len);
			else	putchar('+');
			if(vins[j].pos==p)	putchar(' ');
			if(j+1<vins.size())	for(uint	k=vins[j].pos+1;	k<vins[j+1].pos;	k++)	printf("%s",	k==p?"   ":" ");
		}
		printf("\n");
	}
	else{
		uint	line=0;
		for(uint	i=0;	i<vins.size();	i++)	if(vins[i].str.size()>line)	line=vins[i].str.size();
		for(uint	i=0;	i<line;	i++){
			for(uint	j=s;	j<vins[0].pos;	j++)	printf("%s",	j==p?"   ":" ");
			for(uint	j=0;	j<vins.size();	j++){
				if(vins[j].pos==p)	putchar(' ');
				if(vins[j].len>i)	me->art_char(vins[j].str[i],	'?',	vins[j].qua[i]);
				else	putchar(' ');
				if(vins[j].pos==p)	putchar(' ');
				if(j+1<vins.size())	for(uint	k=vins[j].pos+1;	k<vins[j+1].pos;	k++)	printf("%s",	k==p?"   ":" ");
			}
			printf("\n");
		}
	}		
	return	0;
}

void	SNPView::add_list(cchr	*F){
	FILE	*f=fopen(F,	"rt");
	if(f==NULL){	fprintf(stderr,	"can't open %s!\n",	F);	return;	}
	char	buff[1024],	temp[1024];
	while(fgets(buff,	1024,	f)!=NULL){
		sscanf(buff,	"%s",	temp);	bam_file.push_back(temp);
	}
	fclose(f);
}

void	SNPView::load_ref(cchr	*F){
	gzFile	f=gzopen(F,	"rt");
	if(f==Z_NULL){	fprintf(stderr,	"can't open %s!\n",	F);	return;	}
	char	buff[1024];
	while(gzgets(f,	buff,	1024)!=NULL)	if(buff[0]!='>'){
		for(char	*p=buff;	*p>=33;	p++)	ref.push_back(*p);
	}	
	gzclose(f);
}
	
void	SNPView::document(void){
	fprintf(stderr,	"\nname:	snpview 0.2");
	fprintf(stderr,	"\nfunc:	a text based multiple BAM alignment viewer");
	fprintf(stderr,	"\nauthor:	Yi Wang @ Fudan University");
	fprintf(stderr,	"\nusage:	snpview [options] <chr pos> [bam1 bam2...]");
	fprintf(stderr,	"\n	-l <bams>	add BAMs listed in a file");
	fprintf(stderr,	"\n	-r <ref.fa>	load reference");
	fprintf(stderr,	"\n	-w <width>	half window size");
	fprintf(stderr,	"\n	-R <val>	set max base quality in red (3)");
	fprintf(stderr,	"\n	-Y <val>	set max base quality in yellow (6)");	
	fprintf(stderr,	"\n	-G <val>	set max base quality in green (18)");				
	fprintf(stderr,	"\n	-c		use vt100 color");
	fprintf(stderr,	"\n	-C		collapse insertions");
	fprintf(stderr,	"\n	-v		show variant reads only");
	fprintf(stderr,	"\n	-n		don't reset the console");
	fprintf(stderr,	"\n	-m		filter read by mask");
	fprintf(stderr,	"\n\n");
	exit(0);
}

void	SNPView::options(int	ac,	char	**av){
	wid=32;	red=2;	yel=6;	gre=18;	msk=0;
	isc=false;	col=false;	var=false;	res=true;
		
	char	opt;
	while((opt=getopt(ac,	av,	"l:r:w:cR:Y:G:Cvnm:"))>=0){
		switch(opt) {
		case	'l':	add_list(optarg);	break;
		case	'r':	load_ref(optarg);	break;
		case	'w':	wid=atoi(optarg);	break;
		case	'R':	red=atoi(optarg);	break;
		case	'Y':	yel=atoi(optarg);	break;		
		case	'G':	gre=atoi(optarg);	break;
		case	'c':	isc=true;	break;		
		case	'C':	col=true;	break;						
		case	'v':	var=true;	break;
		case	'n':	res=false;	break;
		case	'm':	msk=atoi(optarg);	break;														
		default:	document();
		}
	}
	if(ac<optind+2)	document();
	//	fill the region for bam_parse_region
	pos=atol(av[optind+1])-1;	min=pos>wid?pos-wid:0;	
	uint	len=ref.size();	max=len?(pos+wid<len?pos+wid:len-1):pos+wid;
	char	temp[256];	
	sprintf(temp,	"%s:%u-%u",	av[optind],	pos+1,	pos+1);
	reg=temp;

	for(int	i=optind+2;	i<ac;	i++)	bam_file.push_back(av[i]);
	sort(bam_file.begin(),	bam_file.end());		
}

void	SNPView::show_bam(cchr	*F){
	printf("%s",	F);	
	uint	len=strlen(F);	
	if(len<2*wid+3){
		for(uint	i=len;	i<2*wid+3;	i++)	putchar(' ');
		printf("\n");
	}
	else	printf("\n");	
	ref_count=0;
	samfile_t	*file=samopen(F, "rb", 0);
	if(file==NULL){	fprintf(stderr,	"can't open %s\n",	F);	return;	}	
	bam_index_t *idx=bam_index_load(F);
	if(idx==NULL){	fprintf(stderr,	"BAM indexing file is not available\n");	return;	}
	int	rid,	beg,	end;
	bam_parse_region(file->header,	reg.c_str(),	&rid,	&beg,	&end);
	if(rid<0){	fprintf(stderr,	"invalid region %s\n",	reg.c_str());	return;	}
	bam_fetch(file->x.bam,	idx,	rid,	beg,	end,	this,	SNPView::fetch);
	bam_index_destroy(idx);	samclose(file);
	if(var&&ref_count)	printf("%u reference reads collapsed\n",	ref_count);
}

int	main(int	ac,	char	**av){
	SNPView	sv;
	sv.options(ac,	av);
	zlibVersion();
	if(sv.res&&system("reset")<0)	printf("\n");
	bool	has_ref=sv.ref.size()>0;
	for(uint	i=sv.min;	i<sv.pos;	i++)	printf("%c",	has_ref?sv.ref[i]:'N');
	printf("|%c|",	has_ref?sv.ref[sv.pos]:'N');
	for(uint	i=sv.pos+1;	i<=sv.max;	i++)	printf("%c",	has_ref?sv.ref[i]:'N');
	printf(" bq mapq\n");
	for(uint	i=0;	i<sv.bam_file.size();	i++)	sv.show_bam(sv.bam_file[i].c_str());	

	return	0;
}

