/*
 ============================================================================
 Name        : 	Encoder.c
 Author      : 	Aníbal Guerra Soler
 Parallel    :  Juan Camilo Peña Vahos
 Copyright   : 	All rights reserved to UdeaCompress
 Description :	Encoder algorithm of the UdeaCompress FASTQ compressor
 ============================================================================
*/

// LIBRERÍAS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <sys/time.h>
#include <omp.h>

// DEFINE
#define NAMES_SIZE 100				// LONGITUD DEL NOMBRE DE LOS ARCHIVOS
#define BASE_BITS 8					// MACROS DEL RADIXSORT
#define BASE (1 << BASE_BITS)
#define MASK (BASE-1)
#define DIGITS(v, shift) (((v) >> shift) & MASK)
#define BYTES_PER_ERROR 2			// 1 BYTE PARA EL OFFSET, 1 BYTE PARA LA DESCRIPCIÓN
#define TEST_PRE	0				// SI TEST_PRE ES 1, SE ACTIVAN LOS ARCHIVOS DE PRUEBA, DE LO CONTRARIO NO
#define ERROR_LOG	0				// SI ERROR_LOG ES 1, SE ACTIVA LA GENERACIÓN DE LOGS DE 
									// ERRORES CONTROLADOS EN LA CODIFICACIÓN
#define TEST_BINST	0				// SI TEST_BINST ES 1, SE ACTIVAN LOS ARCHIVOS DE PRUEBA DE BinINST

// FUNCTIONS PROTOTYPES
//*******************Coding (Compression)(Inst --> binary coding)************************************//
void Inst2Bin(  uint8_t *BinInst, uint8_t *Preambulos, uint32_t *posBInst, uint32_t *posPream,
				char strand, uint8_t MoreFrags, uint16_t lendesc, uint16_t *Offsets, 
				uint8_t *Oper, uint8_t *BaseRead, uint8_t *BaseRef, uint64_t i, uint8_t *flagPream,
				FILE *PREAMBULOS, FILE *ELOGS, FILE *BININST );
uint8_t TrdBitInst( int counter, uint8_t  rest, uint8_t  *Oper, uint8_t  *BaseRead, 
                    uint8_t BaseRef, uint16_t *offset, uint16_t lendesc , char strand, 
                    int *aux_i, FILE *ELOGS);
uint8_t BitsOperR(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, uint16_t lendesc , int *ii);
uint8_t BitsOperF(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, int *ii );
uint8_t Preambulo(uint8_t moreFrags, char strand, uint16_t lendesc, uint8_t flagPream, uint8_t actual);
uint8_t Offset(uint16_t offset, uint8_t *rest);
uint8_t BitsBase(uint8_t BRead, uint8_t BRef, FILE *ELOGS);
void EscalarBases(uint8_t *Base);
void prefix_sum( uint16_t *lendesc, uint32_t *prefixLendesc , uint32_t TotalReads, int NThreads);

//**********************************************SORTING**********************************************//
void RadixSort(uint32_t TotalReads, uint32_t *MapPos, uint64_t *Indexes);

int main(int argc, char *argv[] ) {

	

	// ARGUMENTOS DE ENTRADA
	int 		NThreads = 4;	// NÚMERO DE HILOS POR DEFECTO 4

	// VARIABLES DE PROCESO
	uint32_t	TotalReads;		// CANTIDAD TOTAL DE READS = B*C
	uint64_t	NTErrors;		// CANTIDAD TOTAL DE ERRORES
	uint32_t	B;				// CANTIDAD BASE DE READS
	uint8_t		C;				// COVERAGE DE LA CANTIDAD DE READS
	FILE 		*ALIGN;			// PUNTEROS A LOS ARCHIVOS
	FILE		*PREAMBULOS;	// PUNTERO AL ARCHIVO DE PRUEBA DE PREAMBULOS
	FILE		*ELOGS;			// PUNTERO AL ARCHIVO DE LOGS DE ERRORES
	FILE		*BININST;		// PUNTERO AL ARCHIVO DE PRUEBA DE BININST

	// VARIABLES DE OPERACIÓN
	uint64_t	*Indexes;		// Índices referentes a los Reads
	uint32_t	posBInst;		// Índice que controla BinInst
	uint32_t	*MapPos;        // Posición de Matching respecto a la referencia
	uint16_t  	*lendesc;    	// Cantidad de errores total en el Read
	uint32_t  	*prefixLendesc; // Sumatoria de prefijo exclusivo de lendesc
	char      	*strand;   		// Caractér con el sentido del matching
	uint8_t   	**Oper;			// Arreglo con la operación por error
	uint16_t  	**Offset;   	// Arreglo de offsets por cada error
	uint8_t   	**BaseRef;   	// Arreglo con la base de la referencia (Read Referencia)
	uint8_t 	**BaseRead; 	// Arreglo con la base después de la mutación (Read Destino)

	// VARIABLES DE SALID DEL Inst2Bin
	uint8_t		*BinInst;		// Arreglo de salida del Inst2Bin
	uint8_t		*Preambulos;	// Arreglo de salida con los preámbulos

	// LEER LOS PARÁMETROS DE ENTRADA
	if ( argc > 1 ) {
		for ( int i = 1; i < argc; i++ ) {
			if ( strcmp(argv[i], "-N" ) == 0 ) {
				NThreads	=	(int) atoi (argv[i+1]);
			}
		}
	}

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
		prefixLendesc  =   (uint32_t*)  malloc(TotalReads*sizeof(uint32_t));
		if ( prefixLendesc == NULL ) printf ("Not enough memory for prefixLendesc");
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
		
		// OBTENER TODA LA INFORMACIÓN DESDE EL ARCHIVO DE ALINEAMIENTO
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
		fscanf( ALIGN, "%"SCNu64"", &NTErrors );
		printf("Número de Errores: %"PRIu64"\n",NTErrors);
	}
	fclose (ALIGN);	// SE CIERRA EL ARCHIVO DE ALINEAMIENTO

	

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
	uint64_t TamBinInst	= NTErrors*BYTES_PER_ERROR;
	uint32_t TamPreabulo = floor( TotalReads/2 )+1;
	BinInst	=   (uint8_t*)  malloc(TamBinInst*sizeof(uint8_t));
	if ( BinInst == NULL ) printf ("Not enough memory for BinInst");
	Preambulos	=	(uint8_t*) malloc (TamPreabulo*sizeof(uint8_t));
	if ( Preambulos == NULL ) printf ("Not enough memory for Preambulos");

	posBInst	=	0;
	
	if ( TEST_PRE == 1 || TEST_PRE == 2 ) PREAMBULOS	= fopen( "Preambulos.txt" , "w" );
	if ( ERROR_LOG == 1 ) ELOGS = fopen( "ELogs.txt", "w" );
	if ( TEST_BINST == 1 || TEST_BINST == 2 ) BININST = fopen( "BinInst.txt", "w" );

	printf("NTHREADS %d\n", NThreads);

	// Calcular el arreglo de prefijos exclusivos
	// prefix_sum(lendesc,prefixLendesc,TotalReads, NThreads);

	/*for ( int i = 0; i < TotalReads; i++ ) {
		if( TEST_PRE == 2 ) fprintf(PREAMBULOS, "Lendesc: %"PRIu16" - prefix: %"PRIu32"\n",lendesc[i], prefixLendesc[i]);
	}*/

	// ESTRUCTURA PARA MEDIR TIEMPO DE EJECUCIÓN
	struct timeval t1,t2;
	double elapsedTime;
	gettimeofday(&t1,NULL);

	#pragma omp parallel num_threads(NThreads)
	{
		uint8_t id, Numthreads;
		uint32_t istart, iend, chuncksize;
		Numthreads = omp_get_num_threads();			// NÚMERO DE HILOS CORRIENDO

		chuncksize	=	TotalReads / Numthreads;
		if ( chuncksize % 2 != 0 ) {	// Si es impar
			chuncksize = chuncksize + 1;
		}

		id 		= 	omp_get_thread_num();		// ID DEL HILO
		istart	=	id*chuncksize;				// I INICIAL PARA CADA HILOS
		iend	=	(id+1)*chuncksize;			// I FINAL PARA HILO
		if ( id == Numthreads - 1 ) iend = TotalReads;
		
		printf("Hilo: %d, Start: %d, End: %d\n", id,istart,iend );

		uint32_t posPream;
		if ( id == 0 ) {
			posPream = 0;
		} else {
			posPream = ( (TotalReads/2) / Numthreads ) * id;
		}
		uint8_t flagPream	=	0;
		for ( int index = istart; index < iend; index++ ) {

			// Verificar si el siguiente read mapea en la misma posición
			uint8_t	MoreFrags;
			uint64_t AuxInd	=	Indexes[index];

			// TAREA 1
			if ( (index < TotalReads-1) && (MapPos[AuxInd]	==	MapPos[AuxInd+1]) ) MoreFrags =	1;
			else MoreFrags	=	0;
			
			//Aplicar el inst2bin
			
			Inst2Bin(	BinInst, Preambulos,&posBInst,&posPream, strand[AuxInd],MoreFrags,
						lendesc[AuxInd],Offset[AuxInd],Oper[AuxInd],
						BaseRead[AuxInd],BaseRef[AuxInd],AuxInd, &flagPream, PREAMBULOS, 
						ELOGS, BININST );
			
			// TAREA 6
			if(Offset[AuxInd])		free(Offset[AuxInd]);
			if(Oper[AuxInd]) 		free(Oper[AuxInd]);		
			if(BaseRead[AuxInd]) 	free(BaseRead[AuxInd]);		
			if(BaseRef[AuxInd]) 	free(BaseRef[AuxInd]);	
			
		}
 
	}

	// SE CALCULA EL TIEMPO TOTAL DE EJECUCIÓN Y SE MUESTRA
	gettimeofday(&t2,NULL);
	elapsedTime = (double) (t2.tv_usec - t1.tv_usec) / 1000000 + (double) (t2.tv_sec - t1.tv_sec);
	printf("Processing time: %lf seg\n",elapsedTime);
	printf("Número de Reads: %"PRIu32"\n",TotalReads);
	printf("Número de Errores: %"PRIu64"\n",NTErrors);

	if ( TEST_PRE == 2 ) {
		for ( int i = 0; i < TamPreabulo; i++ ) {
			fprintf(PREAMBULOS,"%"PRIu8"\n", Preambulos[i]);
		}
	}

	if ( TEST_BINST == 2 ) {
		for ( int i = 0; i < TamBinInst; i++ ) {
			fprintf(BININST,"%"PRIu8"\n", BinInst[i]);
		}
	}

	if ( TEST_PRE == 1 || TEST_PRE == 2) fclose(PREAMBULOS);
	if ( ERROR_LOG == 1 ) fclose(ELOGS);
	if ( TEST_BINST == 1 || TEST_BINST == 2) fclose(BININST);

	// CIERRE DE ARCHIVOS Y SE LIBERA LA MEMORIA FALTANTE
	if(MapPos)		free(MapPos);
	if(strand)		free(strand);
	if(lendesc)		free(lendesc);
	if(Offset)		free(Offset);
	if(Oper)		free(Oper);
	if(BaseRead)	free(BaseRead);
	if(BaseRef)		free(BaseRef);
	if(BinInst) 	free(BinInst);
	if(Preambulos) 	free(Preambulos);
	if(Indexes) 	free(Indexes);
    
	

    return 0;

}

// FUNCTIONS

/**
 * @param: BinInst	 -> Arreglo de salida, depende de la cantidad de reads y mutaciones
 * @param: Preambulos-> Arreglo de salida, corresponde a los preámbulos de los reads
 * @param: posBInst  -> Índice de BinInst
 * @param: posPream	 -> Índice de Preambulos
 * @param: strand	 -> Sentido del matching (Forward(F), Reverse(R), Complement(C), Reverse Complement(E))
 * @param: MoreFrags -> Bandera que indica si el siguiente read mapea en la misma posición
 * @param: lendesc	 -> Cantidad de errores del read
 * @param: Offsets	 -> Vector de offsets entre errores
 * @param: Oper		 -> Vector de operaciones
 * @param: BaseRead	 -> Vector de Bases en el Read
 * @param: BaseRef	 -> Vector de Bases en la referencia
 * @param: Index	 -> Posición de este read de acuerdo al nuevo ordenamiento
 * @param: flagPream -> Indica cuando pasar a la siguiente posición en los preámbulos
*/ 
void Inst2Bin(  uint8_t *BinInst, uint8_t *Preambulos, uint32_t *posBInst, uint32_t *posPream, 
				char strand, uint8_t MoreFrags, uint16_t lendesc, uint16_t *Offsets, 
				uint8_t *Oper, uint8_t *BaseRead, uint8_t *BaseRef, uint64_t Index, uint8_t *flagPream,
				FILE *PREAMBULOS, FILE *ELOGS, FILE *BININST){

	uint32_t    auxPosInst =   *posBInst ;
	uint32_t    auxPosPream	=	*posPream;
	uint8_t     rest    =   0x0; 
    uint8_t     aux 	=   0;
    uint8_t     MoreErr =   1;
	int aux_i;
	if ( (*flagPream) == 0 ) {
		// En este caso llena los 4 bits más significativos
		Preambulos[auxPosPream] = 	0;
		Preambulos[auxPosPream]	=	Preambulo(MoreFrags,strand,lendesc,*flagPream,Preambulos[auxPosPream]);
		if ( TEST_PRE == 1 ) fprintf(PREAMBULOS,"Strand: %c , Preambulo: %"PRIu8", Position: %"PRIu32" \n", strand,Preambulos[auxPosPream],auxPosPream);
		(*flagPream) = 1;
	} else {
		// En este caso llena los 4 bits menos significativos
		Preambulos[auxPosPream]	=	Preambulo(MoreFrags,strand,lendesc,*flagPream,Preambulos[auxPosPream]);		
		if ( TEST_PRE == 1 ) fprintf(PREAMBULOS,"Strand: %c , Preambulo: %"PRIu8", Position: %"PRIu32" \n", strand,Preambulos[auxPosPream],auxPosPream);
		(*flagPream) = 0;
		auxPosPream++;
	}


	/*auxPosInst++;
	BinInst[auxPosInst] =   Preambulo(MoreFrags,strand,lendesc);*/
	
    /*if ( lendesc > 0 ){
        if ((strand=='r')||(strand=='e')){
			for (uint8_t  u=0; u<lendesc; u++){ //Converting each separated error of the read
				auxPosInst++;
				BinInst[auxPosInst] = Offset(Offsets[u], &rest);
				if ( TEST_BINST == 1 ) fprintf(BININST,"Offset R (BYTE 1) Index: %"PRIu64" AuxPosInst: %"PRIu32", BinInst[auxPosInst]: %"PRIu8", Offset: %"PRIu16", Rest: %"PRIu8"\n",Index,auxPosInst,BinInst[auxPosInst],Offsets[u],rest);
				auxPosInst++;

				if ( (Oper[u] == 's') || (Oper[u] == 'S') || (Oper[u] == 'i') ) {
					EscalarBases(&BaseRead[u]);
					EscalarBases(&BaseRef[u]);
				}
				BinInst[auxPosInst] = TrdBitInst(u, rest, Oper, BaseRead, BaseRef[u], Offsets, lendesc, strand, &aux_i, ELOGS);
				if ( TEST_BINST == 1 ) fprintf(BININST,"Offset R (BYTE 2) Index: %"PRIu64" AuxPosInst: %"PRIu32", BinInst[auxPosInst]: %"PRIu8", BaseRef: %"PRIu8" \n",Index,auxPosInst,BinInst[auxPosInst],BaseRef[u]);
				
				u=aux_i;

				// if ((BaseRead[u]>=0)&&(BaseRead[u]<=4)){} // FOR REAL ALIGNERS - BOTH CASES					
			}
		}else{         
            for (int  u=lendesc-1;u>=0; u--){				
				auxPosInst++;
				BinInst[auxPosInst]= Offset(Offsets[u+1], &rest);
				if ( TEST_BINST == 1 ) fprintf(BININST,"Offset F (BYTE 1) Index: %"PRIu64" AuxPosInst: %"PRIu32", BinInst[auxPosInst]: %"PRIu8", Offset: %"PRIu16", Rest: %"PRIu8"\n",Index,auxPosInst,BinInst[auxPosInst],Offsets[u],rest);
				auxPosInst++;

				if ( (Oper[u] == 's') || (Oper[u] == 'S') || (Oper[u] == 'i') ) {
					EscalarBases(&BaseRead[u]);
					EscalarBases(&BaseRef[u]);
				}
				BinInst[auxPosInst]= TrdBitInst(u, rest, Oper, BaseRead, BaseRef[u], Offsets, lendesc, strand,&aux_i, ELOGS);
				if ( TEST_BINST == 1 ) fprintf(BININST,"Offset R (BYTE 2) Index: %"PRIu64" AuxPosInst: %"PRIu32", BinInst[auxPosInst]: %"PRIu8", BaseRef: %"PRIu8" \n",Index,auxPosInst,BinInst[auxPosInst],BaseRef[u]);

				u=aux_i;			
			}
		}
    }
	*posBInst=auxPosInst;*/
	*posPream=auxPosPream;
};

uint8_t Preambulo(uint8_t moreFrags, char strand, uint16_t  lendesc, uint8_t flagPream, uint8_t actual){

	uint8_t mask	=	0x01; 
	uint8_t aux		=	actual;

	if (moreFrags==1) {
		aux=mask|aux;   aux=aux<<3;
	}
	//printf ("lendesc %i EDis %i ///", lendesc, EdDis);
	if ((lendesc==0)){ //PERFECT MATCH CUIDADOOOOOOOOOOOOOOOOO
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
    /*if ((lendesc==1)){
			if (strand=='R') mask=0x5;//CreateMask8B(3,1)|CreateMask8B(1,1);//0x101;   //Reverse
			if (strand=='E') mask=0x7; //CreateMask8B(3,3);//0x111;   //rEvErsE complEmEnt
    }*/

	if ( flagPream == 0 ) {
		// 4 bits mas significativos
		aux	=	mask | aux;   
		aux =	aux << 4;
	} else {
		// 4 bits menos significativos
		aux = 	mask|aux;
	}

	
	return(aux);
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
                    char strand, int *aux_i, FILE *ELOGS) {

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
			mask    =   BitsBase(BaseRead[ii], BaseRef, ELOGS);
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
uint8_t BitsBase(uint8_t BRead, uint8_t BRef, FILE *ELOGS){
	
    //calcula la distancia entre la base de la referencia y la base del Read
	//VectorCircular 0:A 1:C 2:G 3:T  4:N
	uint8_t aux =   0x0;
	int auxInt;

	if ((BRead!=0x9)&&(BRef!=0x9)){

		if (BRead>BRef){
			auxInt  =   abs((int)BRead-(int)BRef);
		}else{
			if (BRead==BRef) {
	
				if ( ERROR_LOG == 1 )  fprintf(ELOGS," A Error Grave entre bases iguales Base1 %"PRIu8" Base2 %"PRIu8" \n", BRead,  BRef);
			} else {
				auxInt=((5-(int)BRef)+(int)BRead);
			}
		}
		switch(auxInt){
			case 0:
					
				if ( ERROR_LOG == 1 )  fprintf(ELOGS," B Error Grave entre bases iguales Base1 %"PRIu8" Base2 %"PRIu8" \n", BRead,  BRef);
			break;
			case 1: aux=0x0;  //Distancia 1
			break;
			case 2: aux=0x1;//CreateMask8B(1,1);  //Distancia 2 0x01
			break;
			case 3: aux=0x2;//CreateMask8B(2,1);  //Distancia 3 0x10
			break;
			case 4: aux=0x3;//CreateMask8B(2,2); //Distancia 4 0x11
			break;
			default:
	
				if ( ERROR_LOG == 1 )  fprintf(ELOGS,"Error en la conversión circular Base1 %"PRIu8" Base2 %"PRIu8" \n",  BRead,  BRef);
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
void RadixSort(uint32_t TotalReads, uint32_t *MapPos, uint64_t *Indexes) {

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

/* EscalarBases: Cuando strand sea s S i, se cambian las bases así:
 * A -> 0
 * C -> 1
 * G -> 2
 * T -> 3
 * N -> 4
 */ 
void EscalarBases(uint8_t *Base) {
	switch (*Base){
		case 'A':
			*Base = 0;
		break;
		case 'C':
			*Base = 1;
		break;
		case 'G':
			*Base = 2;
		break;
		case 'T':
			*Base = 3;
		break;
		case 'N':
			*Base = 4;
		break;
	}
}

void prefix_sum( uint16_t *lendesc, uint32_t *prefixLendesc , uint32_t TotalReads, int NThreads) {

	uint32_t nthr, *z;
	uint32_t *x = prefixLendesc;

  	#pragma omp parallel num_threads(NThreads)
  	{
    	int i;
    	#pragma omp single
    	{
      		nthr = omp_get_num_threads();
      		z = malloc(sizeof(uint32_t)*nthr+1);
      		z[0] = 0;
    	}

    	int tid = omp_get_thread_num();
    	uint32_t sum = 0;

    	#pragma omp for schedule(static) 
    		for(i=0; i< TotalReads; i++) {
      			sum += lendesc[i];
      			x[i] = sum;
    		}
    	z[tid+1] = sum;
    	#pragma omp barrier

    	int offset = 0;
    	for(i=0; i<(tid+1); i++) {
        	offset += z[i];
    	}

    	#pragma omp for schedule(static)
    	for(i=0; i< TotalReads; i++) {
      		x[i] += offset;
    	}
  	}
  	free(z);
}
