#include	<stdint.h>
union	XSA_union_double{	uint64_t	i;	double	f;	};
uint64_t	XSA_Seeds[2];
inline	void	srand64(uint64_t	S){
	XSA_Seeds[0]=S*0x5DEECE66Dull+1013904223ull;	XSA_Seeds[1]=0x1234567890ABCDEFull;
}

inline	uint64_t	rand64(void){
	uint64_t	s1=XSA_Seeds[0],	s0=XSA_Seeds[1];	XSA_Seeds[0]=s0;	s1^=s1<<23;	
	return	(XSA_Seeds[1]=(s1^s0^(s1>>17)^(s0>>26)))+s0;
}

double	drand64(void){
	union	XSA_union_double	u;	u.i=(rand64()&0xfffffffffffffull)|0x3ff0000000000000ull;	
	return	u.f-1.0;
}
