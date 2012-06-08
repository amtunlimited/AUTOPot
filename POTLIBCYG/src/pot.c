#include "mathlink.h"
#define N3ATOM 75
#define NATOM 25
#define ISURF 5
#define JSURF 15

extern struct {
	double CART[3][NATOM];
	double ANUZERO;
	int NULBL[NATOM];
	int NFLAG[20];
     int NASURF[ISURF+1][ISURF+1];
     int NDER;
} usricm_;

extern struct {
	double PENGYGS;
	double PENGYES[ISURF];
	double PENGYIJ[JSURF];
	double DGSCART[3][NATOM];
	double DESCART[ISURF][3][NATOM];
	double DIJCART[JSURF][3][NATOM];
} usrocm_;

void pot(
	double cart[], int cartLength,
	double anuzero,
	int nulbl[], int nulblLength,
	int nflag[], int nflagLength,
	int nasurf[], int nasurfLength,
	int nder) {

	double output[10];
	double buffer[3];
	int i,j,k;	 

	i=0;
	for (j=0; j<cartLength/3; j++) {
		for (k=0; k<3; k++) {
			usricm_.CART[k][j]=cart[i];
			i++;
		}
	}
	
	usricm_.ANUZERO=anuzero;
	
	for(i=0; i<nulblLength; i++) {
		usricm_.NULBL[i]=nulbl[i];
	}
	
/*	//Set the option flags*/
/*	usricm_.NFLAG[0]=1;*/
/*	usricm_.NFLAG[1]=1;*/
/*	usricm_.NFLAG[2]=0;*/
/*	usricm_.NFLAG[3]=0;*/
/*	usricm_.NFLAG[4]=0;*/
	
	for(i=0; i<nflagLength; i++) {
		usricm_.NFLAG[i]=nflag[i];
	}
	
/*	//clear out nasurf except the first*/
/*	for(i=0; i<6; i++) {*/
/*		for(j=0;j<6;j++) {*/
/*			usricm_.NASURF[i][j]=0;*/
/*		}*/
/*	}*/
/*	*/
/*	usricm_.NASURF[0][0]=1;*/

	i=0;
	for(j=0; j<6; j++) {
		for(k=0; k<6; k++) {
			usricm_.NASURF[k][j]=nasurf[i];
			i++;
		}
	}
	
	//Set the number of derivitaves to 1 (hopefully)
	usricm_.NDER=nder;
	
	pot_();
	
	output[0]=usrocm_.PENGYGS;
	i=1;
	for (j=0; j<3; j++) {
		for (k=0; k<3; k++) {
			output[i]=usrocm_.DGSCART[k][j];
			i++;
		}
	}
	
	MLPutRealList(stdlink, output, 10);

}
int main(int argc, char* argv[]) {
	prepot_();
	return MLMain(argc, argv);
}
