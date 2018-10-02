/*
 ============================================================================
 Name        : 	Encoder.c
 Author      : 	Aníbal Guerra Soler
 Parallel V  :  Juan Camilo Peña Vahos
 Copyright   : 	All rights reserved to UdeaCompress
 Description :	Encoder algorithm of the UdeaCompress FASTQ compressor
 ============================================================================
 */
/*	
 ============================================================================
		FORMATO DEL ARCHIVO DE ALINEAMIENTO
 ============================================================================
    0. B
    1. C
    3. MapPos
    4. lendesc
    5. strand
    6. 0. oper ->SE REPITEN LENDESC VECES
       1. offset ->SE REPITEN LENDESC VECES
       2. baseRef ->SE REPITEN LENDESC VECES
       3. baseRead ->SE REPITEN LENDESC VECES

*/

// LIBRERÍAS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <inttypes.h>

// DEFINE
#define NAMES_SIZE 40	// LONGITUD DEL NOMBRE DE LOS ARCHIVOS

// FUNCTIONS PROTOTYPES
//**********************************Coding (Compression) ( Inst --> binary coding )
void Inst2Bin(  uint8_t *BinInst, uint32_t *posBInst, char strand, uint8_t MoreFrags, 
                uint16_t lendesc, uint16_t *Offsets, uint8_t *Oper, uint8_t *BaseRead, 
                uint8_t *BaseRef, FILE *outB, FILE *outb2c, FILE *bin, long i, uint16_t EdDis);
uint8_t TrdBitInst( int counter, uint8_t  rest, uint8_t  *Oper, uint8_t  *BaseRead, 
                    uint8_t BaseRef, uint16_t *offset, uint16_t lendesc , char strand, 
                    int *aux_i);
uint8_t BitsOperR(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, uint16_t lendesc , int *ii);
uint8_t BitsOperF(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, int *ii );
uint8_t Preambulo(uint8_t moreFrags, char strand, uint16_t lendesc , uint16_t EdDis);
uint8_t Offset( uint16_t offset, uint8_t *rest);
uint8_t BitsBase(uint8_t BRead, uint8_t BRef);

//**********************************File manipulation
void PrintByte2File(FILE *outB, uint8_t num, int nbits);

int main() {

	// VARIABLES DE PROCESO
	uint32_t	TotalReads;		// CANTIDAD TOTAL DE READS = B*C
	uint32_t	B;				// CANTIDAD BASE DE READS
	uint8_t		C;				// COVERAGE DE LA CANTIDAD DE READS
	char 		*RefAlign;		// NOMBRE ARCHIVO ALIGN
	FILE 		*ALIGN;			// PUNTEROS A LOS ARCHIVOS

	// VARIABLES DE OPERACIÓN
	uint32_t	*MapPos;        // Posición de Matching respecto a la referencia
	uint16_t  	*lendesc;    	// Cantidad de errores total en el Read
	char      	*strand;   		// Caractér con el sentido del matching
	uint8_t   	**Oper;			// Arreglo con la operación por error
	uint16_t  	**Offset;   		// Arreglo de offsets por cada error
	uint8_t   	**BaseRef;   		// Arreglo con la base de la referencia (Read Referencia)
	uint8_t 	**BaseRead; 		// Arreglo con la base después de la mutación (Read Destino)
	

	// ALGORITMO PARA LEER LOS DATOS DE LOS ARCHIVOS
	RefAlign	=	(char*) malloc(NAMES_SIZE*sizeof(char));
	if (RefAlign	==	NULL) printf ("Not enough memory for RefAlign");
	strcpy( RefAlign, "GRCh38.align" );
	
	ALIGN	= 	fopen( RefAlign, "r" );
	if( ALIGN != NULL ) {
		fscanf( ALIGN, "%"SCNu32"",&B );
		fscanf( ALIGN, "%"SCNu8"",&C );
		TotalReads	=	( B*C );

		//DECLARACIÓN DE ARREGLOS
		MapPos        =   (uint32_t*)  malloc(TotalReads*sizeof(uint32_t));
		if ( MapPos == NULL ) printf ("Not enough memory for MapPos");
		lendesc        =   (uint16_t*)  malloc(TotalReads*sizeof(uint16_t));
		if ( lendesc == NULL ) printf ("Not enough memory for lendesc");
		strand        =   (char*)  malloc(TotalReads*sizeof(char));
		if ( strand == NULL ) printf ("Not enough memory for strand");
		// ARREGLOS DE ARREGLOS
		Oper        =   (uint8_t**)  malloc(TotalReads*sizeof(uint8_t*));
		if ( Oper == NULL ) printf ("Not enough memory for Oper");
		Offset        =   (uint16_t**)  malloc(TotalReads*sizeof(uint16_t*));
		if ( Offset == NULL ) printf ("Not enough memory for lendesc");
		BaseRef        =   (uint8_t**)  malloc(TotalReads*sizeof(uint8_t*));
		if ( BaseRef == NULL ) printf ("Not enough memory for BaseRef");
		BaseRead        =   (uint8_t**)  malloc(TotalReads*sizeof(uint8_t*));
		if ( BaseRead == NULL ) printf ("Not enough memory for BaseRead");

		printf("%"PRIu32"\n",TotalReads);
		
		for ( int i = 0; i < TotalReads; i++ ) {

			fscanf( ALIGN, "%"SCNu32"", &MapPos[i] );
			fscanf( ALIGN, "%"SCNu16"", &lendesc[i] );
			fscanf( ALIGN, " %c", &strand[i] );

			//printf( "%"PRIu32"\n", MapPos[i] ); fflush(stdout);
			//printf( "%"PRIu16"\n", lendesc[i] ); fflush(stdout);
			//printf( "%c\n", strand[i] ); fflush(stdout);

			if ( lendesc[i] != 0 ) {
				for ( int j = 0; j < lendesc[i]; j++ ) {

					Oper[i]        =   (uint8_t*)  malloc(lendesc[i]*sizeof(uint8_t));
					if ( Oper[i] == NULL ) printf ("Not enough memory for Oper");
					Offset[i]        =   (uint16_t*)  malloc(lendesc[i]*sizeof(uint16_t));
					if ( Offset[i] == NULL ) printf ("Not enough memory for lendesc");
					BaseRef[i]        =   (uint8_t*)  malloc(lendesc[i]*sizeof(uint8_t));
					if ( BaseRef[i] == NULL ) printf ("Not enough memory for BaseRef");
					BaseRead[i]        =   (uint8_t*)  malloc(lendesc[i]*sizeof(uint8_t));
					if ( BaseRead[i] == NULL ) printf ("Not enough memory for BaseRead");

					fscanf( ALIGN, " %c", &Oper[i][j] );
					fscanf( ALIGN, "%"SCNu16"", &Offset[i][j] );
					fscanf( ALIGN, " %c", &BaseRef[i][j] );
					fscanf( ALIGN, " %c", &BaseRead[i][j] );

					//printf( "%c\n", Oper[i][j] );
					//printf( "%"PRIu16"\n", Offset[i][j] );
					//printf( "%c\n", BaseRef[i][j] );
					//printf( "%c\n", BaseRead[i][j] );
				
				}
			}
		}
	}

	fclose( ALIGN );
	
	if(RefAlign)	 free(RefAlign);
    return 0;

}

// FUNCTIONS

/**
 * @param: BinInst	 -> Arreglo de salida, depende de la cantidad de reads y mutaciones
 * @param: posBInst  -> Índice de BinInst
 * @param: strand	 -> Sentido del matching (Forward(F), Reverse(R), Complement(C), Reverse Complement(E))
 * @param: MoreFrags -> Bandera que indica si el siguiente read mapea en la misma posición
 * @param: lendesc	 -> Cantidad de errores del read
 * @param: Offsets	 -> Vector de offsets entre errores
 * @param: Oper		 -> Vector de operaciones
 * @param: BaseRead	 -> Vector de Bases en el Read
 * @param: BaseRef	 -> Vector de Bases en la referencia
 * @param: Index	 -> Posición de este read de acuerdo al nuevo ordenamiento
*/ 
void Inst2Bin(  uint8_t *BinInst,  uint32_t *posBInst, char strand, uint8_t MoreFrags, 
                uint16_t lendesc, uint16_t *Offsets, uint8_t *Oper, uint8_t *BaseRead, 
                uint8_t *BaseRef, FILE *outB, FILE *outb2c, FILE *bin, long Index, uint16_t EdDis){
	
	uint32_t    auxPosInst =   *posBInst ;
	uint8_t     rest    =   0x0; 
    uint8_t     aux =   0;
    uint8_t     MoreErr =   1;

	int aux_i;

	auxPosInst++;
	BinInst[auxPosInst] =   Preambulo(MoreFrags,strand,lendesc,EdDis);

	PrintByte2File(outB, BinInst[auxPosInst], 8); // Print to a Bytes File (as a string of '0' and '1' )

    if ( ((lendesc>0)&&((strand=='r')||(strand=='e'))) || ((lendesc>1)&&((strand=='f')||(strand=='c'))) ){
		
        if ((strand=='r')||(strand=='e')){
			for (uint8_t  u=0; u<lendesc; u++){ //Converting each separated error of the read
				if ((Oper[u]!='_')&&(BaseRead[u]>=0)&&(BaseRead[u]<=4)){

					auxPosInst++;

					BinInst[auxPosInst] = Offset(Offsets[u], &rest);
					PrintByte2File(outB, BinInst[auxPosInst], 8);

					auxPosInst++;

					BinInst[auxPosInst] = TrdBitInst(u, rest, Oper, BaseRead, BaseRef[u], Offsets, lendesc, strand, &aux_i);
					PrintByte2File(outB, BinInst[auxPosInst], 8);

					u=aux_i;
				}
			}
		}else{  
            
            for (int  u=lendesc-1;u>=0; u--){

				if ((Oper[u]!='_')&&(BaseRead[u]>=0)&&(BaseRead[u]<=4)){
					auxPosInst++;

					BinInst[auxPosInst]= Offset(Offsets[u+1], &rest);
					PrintByte2File(outB, BinInst[auxPosInst], 8);

					auxPosInst++;

                    BinInst[auxPosInst]= TrdBitInst(u, rest, Oper, BaseRead, BaseRef[u], Offsets, lendesc, strand,&aux_i);
					PrintByte2File(outB, BinInst[auxPosInst], 8);

					u=aux_i;
				}
			}
		}
    }
	*posBInst=auxPosInst;
	fprintf(outB,"\n \n");
};

uint8_t Preambulo(uint8_t moreFrags, char strand, uint16_t  lendesc, uint16_t EdDis ){

	uint8_t mask    =   0x01;
    uint8_t aux     =   0x0;

	if (moreFrags==1) {
		aux=mask|aux;   aux=aux<<3;
	}
	if ((lendesc<=1)){ 
		if (strand=='F') mask=0x0;  //Forward
		if (strand=='R') mask=0x1;  //Reverse
		if (strand=='C') mask=0x2;  //Complement 10
		if (strand=='E') mask=0x3;  //11 CreateMask8B(2,2);//0x011;   //rEvErsE complEmEnt
	}else{
		if (strand=='F') mask=0x4;  //CreateMask8B(3,1);//0x100;   //Forward
		if (strand=='R') mask=0x5;  //CreateMask8B(3,1)|CreateMask8B(1,1);//0x101;   //Reverse
		if (strand=='C') mask=0x6;  //CreateMask8B(3,2);//0x110;   //Complement
		if (strand=='E') mask=0x7;  //CreateMask8B(3,3);//0x111;   //rEvErsE complEmEnt
	};

    if ((lendesc==1)&&(EdDis==1)){
		if (strand=='R') mask=0x5;  //CreateMask8B(3,1)|CreateMask8B(1,1);//0x101;   //Reverse
		if (strand=='E') mask=0x7;  //CreateMask8B(3,3);//0x111;   //rEvErsE complEmEnt
    }
	aux=mask|aux;
	return(aux);
};

void PrintByte2File(FILE *outB, uint8_t num, int nbits){
	for( int i = nbits-1 ; i>=0 ; --i )
	    if( num & (1 << i) ) fprintf(outB,"%c", '1');
	        else fprintf(outB,"%c", '0');
	 fprintf(outB,"\t");
};

/*  8 bits menos significativos del offset. Los 2 bits + significativos en el 3er Entero, 
    en las posiciones menos significativas*/
uint8_t Offset( uint16_t offset, uint8_t *rest) {

	uint8_t mask=0x01, aux=0x0, auxRest=0x0, mod, bitCnt=0;

    while (offset>0){

    	if (bitCnt<8){
    		mod=offset%2; offset=offset/2;
	    	if (mod==1) aux=aux|mask;
	    	mask=mask<<1;
	    	bitCnt++;
	    	if (bitCnt==8) mask=0x01;
    	}else{
    		mod=offset%2; offset=offset/2;
    		if (mod==1) auxRest=auxRest|mask;
    		mask=mask<<1;
    	}
    }

    if (offset==1)
    	if (bitCnt>=8){
    		auxRest=auxRest|mask;
    	}else{
    		aux=aux|mask;
    	}

    *rest   =   auxRest;
	return(aux);
};

// 2 MSB (aka Most significant bit ): Offset's 2 MSB bits when this is >=256
// 3er bit mas significativo: Bit del MoreError
// 4to, 5to y 6to bit mas significativo: OPER
// 2 bits menos significativos: Base
uint8_t TrdBitInst( int i, uint8_t  rest, uint8_t  *Oper, uint8_t  *BaseRead, 
                    uint8_t BaseRef, uint16_t *offset, uint16_t lendesc, 
                    char strand, int *aux_i) {

	uint8_t mask=0x01, aux=0x0, mask2=0x0;
	int ii  =   i;
	aux =   aux|rest;
    aux =   aux<<1;

    if ((strand=='R')||(strand=='e')) mask2 =   BitsOperR(Oper, BaseRead, offset, lendesc, &i);
    if ((strand=='F')||(strand=='c')) mask2 =   BitsOperF(Oper, BaseRead, offset,  &i);

   	if (((strand=='R')||(strand=='e'))&&((i<lendesc-1)&&(Oper[i+1]!='_'))) aux=mask|aux;  
   		else if (((strand=='F')||(strand=='c'))&&(i>0)) aux=mask|aux;

   	aux =   aux<<3;
    aux =   aux|mask2;
	aux =   aux<<2;

	if((mask2==0x2)||(mask2==0x3)||(mask2==0x4)||(mask2==0x5)||(mask2==0x7)) {
		mask    =   0x00;
	}else {
		if ((mask2  !=  0x1)){
			mask    =   BitsBase(BaseRead[ii], BaseRef);
		}else{
			mask=BaseRead[ii];
		}

		aux=aux|mask;
	}

	*aux_i=i;
	return(aux);
};

uint8_t BitsOperF(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, int *ii ){ 

	uint8_t aux =   0x0, NDelC=1;
	int i   =   *ii;

	switch (oper[i]){

		case 'd': 

			while (((i-1)>=0)&&(offset[i]==0)&&(oper[i-1]=='d')&&(NDelC<=4)){
				NDelC++; 				
                i--;
			}

			if (NDelC==5) { NDelC--;    i++; }

			switch(NDelC){
				case 1:aux=0x2;break;//aux=CreateMask8B(2,1);   // Delecion Simple 0x010;
				case 2:aux=0x3;break;//aux=CreateMask8B(2,2);   //0x011
				case 3:aux=0x4;break;//CreateMask8B(3,1);       //0x100
				case 4:aux=0x5;break;//CreateMask8B(3,1)|CreateMask8B(1,1);  //Delecion Cuadruple Consecutiva //aux=0x101
			}

		break;
		case 'i':

			if( (0x4==baseRead[i])||('n'==baseRead[i]) ) { aux=0x7; }
            else{ aux=0x1; }
			break;

		case 's':

			if(((i-1)>=0)&&(offset[i]==0)&&(oper[i-1]=='s')&&(baseRead[i-1]==baseRead[i])){
                aux=0x6;    //CreateMask8B(3,2);    //=0x110
			    i--;
			}else{ 
				aux=0x0;
			}
		break;
	}
	*ii =   i;
	return (aux);

};


uint8_t BitsOperR(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, uint16_t lendesc , int *ii ){ 

	uint8_t aux=0x0, NDelC=1, i=*ii;

	switch (oper[i]){
		case 'd':

			while ((i<lendesc-1)&&(offset[i+1]==0)&&(oper[i+1]=='d')&&(NDelC<=4)){
				NDelC++;    i++;
			}
			if (NDelC==5) { NDelC--;    i--; }
			switch(NDelC){
				case 1:aux=0x2; break; // Delecion Simple 0x010; CreateMask8B(2,1)
				case 2:aux=0x3; break;//0x011 CreateMask8B(2,2)
				case 3:aux=0x4; break;//0x100 CreateMask8B(3,1)
				case 4:aux=0x5; break; //Delecion Cuadruple Consecutiva //101 CreateMask8B(3,1)|CreateMask8B(1,1)
			}
		break;
		case 'i':
			if(0x4==baseRead[i]) { aux=0x7; } 
            else { aux=0x1 ; }
		break;
		case 's':
			if((i<lendesc-1)&&(offset[i+1]==0)&&(oper[i+1]=='s')&&(baseRead[i]==baseRead[i+1])){
				aux=0x6;
				i++;
			} else { aux=0x0; }
		break;
	}
	*ii =   i;
	return (aux);
};

//Distancia De la BASE en la REF a la Base en el READ, devuelve un numero del 0 al 3,
//al cual eventualmente tocara sumarle 1
uint8_t BitsBase(uint8_t BRead, uint8_t BRef){
	
    //calcula la distancia entre la base de la referencia y la base del Read
	//VectorCircular 0:A 1:C 2:G 3:T  4:N
	uint8_t aux =   0x0;
	int auxInt;

	if ((BRead!=0x9)&&(BRef!=0x9)){

		if (BRead>BRef){
			auxInt  =   abs((int)BRead-(int)BRef);
		}else{
			if (BRead==BRef) printf("Error Grave entre bases iguales Base1 %u Base2 %u \n", BRead,  BRef);
			else{
				auxInt=((5-(int)BRef)+(int)BRead);
			}
		}
		switch(auxInt){
			case 0: printf("Error Grave entre bases iguales Base1 %u Base2 %u \n", BRead,  BRef);
			break;
			case 1: aux=0x0;  //Distancia 1
			break;
			case 2: aux=0x1;//CreateMask8B(1,1);  //Distancia 2 0x01
			break;
			case 3: aux=0x2;//CreateMask8B(2,1);  //Distancia 3 0x10
			break;
			case 4: aux=0x3;//CreateMask8B(2,2); //Distancia 4 0x11
			break;
			default: printf("Error Diferencia Grande bases Base1 %u BAse2 %u \n",  BRead,  BRef);
		}
	}
	return(aux);
};

