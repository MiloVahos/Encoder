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
#include <omp.h>

// DEFINE
#define NAMES_SIZE 40				// LONGITUD DEL NOMBRE DE LOS ARCHIVOS
#define BASE_BITS 8					// MACROS DEL RADIXSORT
#define BASE (1 << BASE_BITS)
#define MASK (BASE-1)
#define DIGITS(v, shift) (((v) >> shift) & MASK)
#define BYTES_PER_ERROR 2			// 1 BYTE PARA EL OFFSET, 1 BYTE PARA LA DESCRIPCIÓN

// FUNCTIONS PROTOTYPES
//**********************************Coding (Compression)(Inst --> binary coding)*********//
void Inst2Bin(  uint8_t *BinInst, uint32_t *posBInst, char strand, uint8_t MoreFrags, 
                uint16_t lendesc, uint16_t *Offsets, uint8_t *Oper, uint8_t *BaseRead, 
                uint8_t *BaseRef, uint64_t i);
uint8_t TrdBitInst( int counter, uint8_t  rest, uint8_t  *Oper, uint8_t  *BaseRead, 
                    uint8_t BaseRef, uint16_t *offset, uint16_t lendesc , char strand, 
                    int *aux_i);
uint8_t BitsOperR(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, uint16_t lendesc , int *ii);
uint8_t BitsOperF(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, int *ii );
uint8_t Preambulo(uint8_t moreFrags, char strand, uint16_t lendesc);
uint8_t Offset(uint16_t offset, uint8_t *rest);
uint8_t BitsBase(uint8_t BRead, uint8_t BRef);

//**********************************File manipulation************************************//
//void PrintByte2File(FILE *outB, uint8_t num, int nbits);

//**********************************SORTING**********************************************//
void RadixSort(int32_t TotalReads, uint32_t *MapPos, uint64_t *Indexes);

int main() {

	// VARIABLES DE PROCESO
	uint32_t	TotalReads;		// CANTIDAD TOTAL DE READS = B*C
	uint64_t	NTErrors;		// CANTIDAD TOTAL DE ERRORES
	uint32_t	B;				// CANTIDAD BASE DE READS
	uint8_t		C;				// COVERAGE DE LA CANTIDAD DE READS
	FILE 		*ALIGN;			// PUNTEROS A LOS ARCHIVOS

	// VARIABLES DE OPERACIÓN
	uint64_t	*Indexes;		// Índices referentes a los Reads
	uint32_t	posBInst;		// Índice que controla BinInst
	uint8_t		MoreFrags;		// Indica si el siguiente Read Mapea en la misma posición
	uint32_t	*MapPos;        // Posición de Matching respecto a la referencia
	uint16_t  	*lendesc;    	// Cantidad de errores total en el Read
	char      	*strand;   		// Caractér con el sentido del matching
	uint8_t   	**Oper;			// Arreglo con la operación por error
	uint16_t  	**Offset;   	// Arreglo de offsets por cada error
	uint8_t   	**BaseRef;   	// Arreglo con la base de la referencia (Read Referencia)
	uint8_t 	**BaseRead; 	// Arreglo con la base después de la mutación (Read Destino)

	// VARIABLES DE SALID DEL Inst2Bin
	uint8_t		*BinInst;		// Arreglo de salida del Inst2Bin

	// 1. OBTENER LOS DATOS QUE PROVIENEN DEL ARG
	ALIGN	= 	fopen( "GRCh38.align" , "r" );
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
		
		for ( int i = 0; i < TotalReads; i++ ) {

			fscanf( ALIGN, "%"SCNu32"", &MapPos[i] );
			fscanf( ALIGN, "%"SCNu16"", &lendesc[i] );
			fscanf( ALIGN, " %c", &strand[i] );

			if ( lendesc[i] != 0 ) {

				Oper[i]        =   (uint8_t*)  malloc(lendesc[i]*sizeof(uint8_t));
				if ( Oper[i] == NULL ) printf ("Not enough memory for Oper");
				Offset[i]        =   (uint16_t*)  malloc(lendesc[i]*sizeof(uint16_t));
				if ( Offset[i] == NULL ) printf ("Not enough memory for lendesc");
				BaseRef[i]        =   (uint8_t*)  malloc(lendesc[i]*sizeof(uint8_t));
				if ( BaseRef[i] == NULL ) printf ("Not enough memory for BaseRef");
				BaseRead[i]        =   (uint8_t*)  malloc(lendesc[i]*sizeof(uint8_t));
				if ( BaseRead[i] == NULL ) printf ("Not enough memory for BaseRead");

				for ( int j = 0; j < lendesc[i]; j++ ) {

					uint8_t oper, baseref, baseread;
					uint16_t offset;

					fscanf( ALIGN, " %c", &oper );
					fscanf( ALIGN, "%"SCNu16"", &offset );
					fscanf( ALIGN, " %c", &baseref );
					fscanf( ALIGN, " %c", &baseread );

					memcpy( &Oper[i][j], &oper, sizeof(uint8_t));
					memcpy( &Offset[i][j], &offset, sizeof(uint16_t));
					memcpy( &BaseRef[i][j], &baseref, sizeof(uint8_t));
					memcpy( &BaseRead[i][j], &baseread, sizeof(uint8_t));
				}
			} else {
				Oper[i]        =   (uint8_t*)  malloc(sizeof(uint8_t));
				if ( Oper[i] == NULL ) printf ("Not enough memory for Oper");
				Offset[i]        =   (uint16_t*)  malloc(sizeof(uint16_t));
				if ( Offset[i] == NULL ) printf ("Not enough memory for lendesc");
				BaseRef[i]        =   (uint8_t*)  malloc(sizeof(uint8_t));
				if ( BaseRef[i] == NULL ) printf ("Not enough memory for BaseRef");
				BaseRead[i]        =   (uint8_t*)  malloc(sizeof(uint8_t));
				if ( BaseRead[i] == NULL ) printf ("Not enough memory for BaseRead");
			}
		}
		fscanf( ALIGN, "%"SCNu64"",&NTErrors );
	}

	// 2. USANDO EL RADIX SORT SE ORDENA EL VECTOR DE ÍNDICES DE ACUERDO CON LA POSICIÓN DE MAPEO
	// 		- SE CREA EL VECTOR DE ÍNDICES [0 - TotalReads-1]
	Indexes	=   (uint64_t*)  malloc(TotalReads*sizeof(uint64_t));
	if ( Indexes == NULL ) printf ("Not enough memory for Indexes");
	for ( int i = 0; i < TotalReads; i++ ) Indexes[i] =	i;

	//		- ALGORITMO DE ORDENAMIENTO RADIX SORT
	uint32_t *AuxMapPos;
	AuxMapPos	=	(uint32_t*) malloc( TotalReads*sizeof(uint32_t));
	memcpy(AuxMapPos,MapPos,TotalReads*sizeof(uint32_t));
	RadixSort(TotalReads,AuxMapPos,Indexes);
	free(AuxMapPos);
	
	//		- APLICACIÓN DEL INS2BIN
	BinInst	=   (uint8_t*)  malloc((TotalReads*NTErrors*BYTES_PER_ERROR)*sizeof(uint8_t));
	if ( BinInst == NULL ) printf ("Not enough memory for Indexes");
	posBInst	=	0;
	MoreFrags	=	0;
	uint64_t AuxInd	=	0;

	for ( int index = 0; index < 10; index++ ) {

		// Verificar si el siguiente read mapea en la misma posición
		if ( MapPos[Indexes[index]]	==	MapPos[Indexes[index+1]] ) MoreFrags	=	1;
		else MoreFrags	=	0;
		
		//Aplicar el inst2bin
		AuxInd	=	Indexes[index];
		Inst2Bin(	BinInst,&posBInst,strand[AuxInd],MoreFrags,
					lendesc[AuxInd],Offset[AuxInd],Oper[AuxInd],
					BaseRead[AuxInd],BaseRef[AuxInd],AuxInd );
		
	}

	//for ( int i = 0; i < TotalReads*NTErrors*BYTES_PER_ERROR; i++ ) printf("%"PRIu8" ",BinInst[i]);
	

	fclose( ALIGN );
	
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
                uint8_t *BaseRef, uint64_t Index){
	uint32_t    auxPosInst =   *posBInst ;
	uint8_t     rest    =   0x0; 
    uint8_t     aux =   0;
    uint8_t     MoreErr =   1;

	int aux_i;

	auxPosInst++;
	BinInst[auxPosInst] =   Preambulo(MoreFrags,strand,lendesc);
    if ( ((lendesc>0)&&((strand=='r')||(strand=='e'))) || ((lendesc>1)&&((strand=='f')||(strand=='c'))) ){
		
        if ((strand=='r')||(strand=='e')){
			for (uint8_t  u=0; u<lendesc; u++){ //Converting each separated error of the read
				if ((Oper[u]!='_')&&(BaseRead[u]>=0)&&(BaseRead[u]<=4)){

					auxPosInst++;

					BinInst[auxPosInst] = Offset(Offsets[u], &rest);

					auxPosInst++;

					BinInst[auxPosInst] = TrdBitInst(u, rest, Oper, BaseRead, BaseRef[u], Offsets, lendesc, strand, &aux_i);

					u=aux_i;
				}
			}
		}else{  
            
            for (int  u=lendesc-1;u>=0; u--){

				if ((Oper[u]!='_')&&(BaseRead[u]>=0)&&(BaseRead[u]<=4)){
					auxPosInst++;

					BinInst[auxPosInst]= Offset(Offsets[u+1], &rest);

					auxPosInst++;

                    BinInst[auxPosInst]= TrdBitInst(u, rest, Oper, BaseRead, BaseRef[u], Offsets, lendesc, strand,&aux_i);

					u=aux_i;
				}
			}
		}
    }
	*posBInst=auxPosInst;
};

uint8_t Preambulo(uint8_t moreFrags, char strand, uint16_t  lendesc){

	uint8_t mask	=	0x01; 
	uint8_t aux		=	0x0;

	if (moreFrags==1) {
		aux=mask|aux;   aux=aux<<3;
	}
	//printf ("lendesc %i EDis %i ///", lendesc, EdDis);
	if ((lendesc<=1)){ //PERFECT MATCH CUIDADOOOOOOOOOOOOOOOOO
		if (strand=='F') mask=0x0; //Forward
		if (strand=='R') mask=0x1;   //Reverse
		if (strand=='C') mask=0x2;   //Complement 10
		if (strand=='E') mask=0x3;//11 CreateMask8B(2,2);//0x011;   //rEvErsE complEmEnt
		//if (strand=='T') mask=0x011; //another Transform
	}else{ //ERROR MATCH
		if (strand=='f') mask=0x4;//CreateMask8B(3,1);//0x100;   //Forward
		if (strand=='r') mask=0x5;//CreateMask8B(3,1)|CreateMask8B(1,1);//0x101;   //Reverse
		if (strand=='c') mask=0x6;//CreateMask8B(3,2);//0x110;   //Complement
		if (strand=='e') mask=0x7; //CreateMask8B(3,3);//0x111;   //rEvErsE complEmEnt
		//if (strand=='T') mask=0x011; //another Transform
	};
    if ((lendesc==1)){
			if (strand=='R') mask=0x5;//CreateMask8B(3,1)|CreateMask8B(1,1);//0x101;   //Reverse
			if (strand=='E') mask=0x7; //CreateMask8B(3,3);//0x111;   //rEvErsE complEmEnt
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
			printf("%"PRIu8",%"PRIu8",%"PRIu8"   ",Oper[i],BaseRead[ii], BaseRef);
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

uint8_t BitsBase(uint8_t BRead, uint8_t BRef){
	
    //calcula la distancia entre la base de la referencia y la base del Read
	//VectorCircular 0:A 1:C 2:G 3:T  4:N
	uint8_t aux =   0x0;
	int auxInt;

	if ((BRead!=0x9)&&(BRef!=0x9)){

		if (BRead>BRef){
			auxInt  =   abs((int)BRead-(int)BRef);
		}else{
			if (BRead==BRef) printf(" A Error Grave entre bases iguales Base1 %u Base2 %u \n", BRead,  BRef);
			else{
				auxInt=((5-(int)BRef)+(int)BRead);
			}
		}
		switch(auxInt){
			case 0: printf(" B Error Grave entre bases iguales Base1 %u Base2 %u \n", BRead,  BRef);
			break;
			case 1: aux=0x0;  //Distancia 1
			break;
			case 2: aux=0x1;//CreateMask8B(1,1);  //Distancia 2 0x01
			break;
			case 3: aux=0x2;//CreateMask8B(2,1);  //Distancia 3 0x10
			break;
			case 4: aux=0x3;//CreateMask8B(2,2); //Distancia 4 0x11
			break;
			default: printf("Error en la conversión circular Base1 %u BAse2 %u \n",  BRead,  BRef);
		}
	}
	return(aux);
};

/*	RadixSort:	Algoritmo de ordenamiento por cubetas, se encargará de ordenar los índices
 *  			de acuerdo con la posición de mapeo del read, MapPos
 * @param:	TotalReads	->	Cantidad de reads que va a procesar el algoritmo
 * @param:	MapPos		->	Vector que contiene las posiciones de mapeo de cada read
 * @param:	Indexes		->	Vector de índices que se van a ordenar de acuerdo a MapPos
*/				
void RadixSort(int32_t TotalReads, uint32_t *MapPos, uint64_t *Indexes) {

	// Buffer de ordenamiento temporal para MapPos
    uint64_t *buffer = (uint64_t *) malloc(TotalReads*sizeof(uint64_t));
    if (buffer  == NULL) printf("No hay espacio suficiente para buffer\n");

	// Buffer de ordenamiento temporal para los índices
	uint32_t *bufferAux = (uint32_t *) malloc(TotalReads*sizeof(uint32_t));
    if (bufferAux  == NULL) printf("No hay espacio suficiente para bufferAux\n");
    
    int total_digits = sizeof(uint64_t)*8;
    int32_t i;
    int shift, cur_t;

    for( shift = 0;	shift < total_digits;	shift += BASE_BITS ) { 

        int64_t bucket[BASE] = {0};
        int64_t local_bucket[BASE] = {0};
        
		for(i = 0; i < TotalReads; i++) local_bucket[DIGITS(MapPos[i], shift)]++;
        for(i = 0; i < BASE; i++) bucket[i] += local_bucket[i];
        for(i = 1; i < BASE; i++) bucket[i] += bucket[i - 1];
            
        int nthreads = 1 ;
        int tid = 0;  
        for(cur_t = nthreads - 1; cur_t >= 0; cur_t--) {
            if(cur_t == tid) {
                for(i = 0; i < BASE; i++) {
                    bucket[i] -= local_bucket[i];
                    local_bucket[i] = bucket[i];
                }
            }
        }
        for(i = 0; i < TotalReads; i++){
			int Index = local_bucket[DIGITS(MapPos[i], shift)]++;
			buffer[Index] = MapPos[i];
			bufferAux[Index] = Indexes[i];	
		} 

		uint32_t* tmp = Indexes;
		int32_t h=0;
		for( i = 0; i < TotalReads; i++ ){
			MapPos[i] = buffer[h];
			Indexes[i] = bufferAux[h]; 
			h++;
		}
	
		h=0;
		for( i = 0; i < TotalReads; i++ ){
			buffer[h] = tmp[i];
			h++;
		}
    }
    if(buffer) free(buffer);
}

