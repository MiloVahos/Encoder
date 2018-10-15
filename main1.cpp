
// 	This file implements the interfaces for the modules that are progressively developed
//	It does not necessarily corresponds to the final user interface

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <fstream>
#include <inttypes.h>
#include "lib/libbzip2/bzlib.h"
#include "lib/MINE/bammatch.h"
#include "lib/SAMTOOLS/sam.h"

//Constants to define/adjust
#define BZIP2_BUFFER_MAX_LENGTH 100
#define UNCOMPRESSED_MAX_LENGTH 100

#ifndef INSTR_LENGTH			//Maximum Length of the (intermediate) Instruction representation built using the alignment info.
	#define INSTR_LENGTH 511 	//All info of each error (op, base, etc) is placed in a single position.
#endif							//This should be related to Read lentgth since the number of errors is a % oh the read length (Max 40%)

#ifndef NUM_READSM				//Number of reads (mapped or not), this value should come from the aligner
	#define NUM_READSM 10000  	//***This value must be higher, (until 10-20.000K)
#endif

#ifndef MAX_ERROR_READ     		//Maximum number of errors per read,
	#define MAX_ERROR_READ 300  //This value must be defined as a percentage of the read length (Max 40%). Must be equal to INSTR_LENGTH
#endif

#ifndef BYTES_PER_ERROR     	//Maximum number of of bytes required to store an error
	#define BYTES_PER_ERROR 3   //This value is OK
#endif

#ifndef REF_LENGTH				//Maximum Length of the reference, in this case a single reference is consider.
								//But in the future, the reference could be a group of chromosomes (in a single or even in separated files)
	#define REF_LENGTH 50000  	//This value should be Much More, around 10.000 K
#endif							//Please set CHRM_LENGTH to the same value in BAMMATCH

#ifndef READ_LENGTH				//Length of any of the reads in the fastqFile.
	#define READ_LENGTH 1024	//This alue is OK. Maximum supported: 2048
#endif							//Please set LINE_LENGTH to the same value in BAMMATCH

#ifndef DEBUG					//Set Debug Mode On=1/Off=0
	#define DEBUG 1
#endif



//#ifndef STR_LENTGH            //YA SE ME OLVIDO!
//#define STR_LENTGH 100
//#endif

//#ifndef INSTR_LENGTH
//	#define INSTR_LENGTH 1024
//#endif


struct Instr_data{ 		//Object structure for the instruction data. This could disappear if the POO focus is discarded
	char type;	 	//'s'ubstitution, 'd'eletion, 'i'nsertion
	char ref_ltr;//Letter placed in the reference position
	char rd_ltr; //Letter placed in the read position
	int repeat;  //This field has a double purpose, 1st: all vallid instructions should have at least a value of 1 here.
				 //If this value is exactly 1 it is a valid instruction. If it is higher, it means that there are several
				 //consecutive 'N' (for ins and and subst, or consecutive deletions)
	int offset;  //El offset.
};

typedef struct Instr_data Instr;

struct Oligo_info{ 		//Object structure for the alignment data.
	char *prefix;       //*Non-used, comes from aligner
	char *suffix;		//*Non-used, comes from aligner
	char *read;		  //Read Sequence
	char *qvals;	  //Quality scores
	uint32_t loc;		  //Position
	char strand; 	  //R/F
	char *descriptor; //MD String
	int EditDist; //Length of the descriptor, calculated by this program (not the aligner). It is the exact value
	int lendesc; 	  ////Length of the descriptor, calculated by this program (not the aligner) when parsing
	int lonRead; 	  //Read length according to the descriptor
	int NREAD; 	  //Number of this read (Work as an Id also)
	int NDel;       //Number of Del operations for this read
	int NIns;       //Number of Ins operations for this read
};

typedef struct Oligo_info Olg;

struct rec  //Aux structure for data type conversion to store in the final bin file
{
	uint8_t x;
};

/****Protypes****/

//**********************************PARSING SAM
Instr *parse_line(Olg olg, int *len, char *ref); 	//SobreElDescriptor
int VerificRead(Olg olg, Instr *MInst, char *ref1, int lenR, uint8_t Match, uint16_t *Of, unsigned char *Op, uint8_t *Ba); //Reconstruye y verifica igualdad del read Mapeado
void add_new_instruction(Instr *instr_ar, int i, char type, char ref_ltr, char rd_ltr,int repeat, int offset);
Olg MIparse_input(bam1_t *b, char *ref, int ref_len, long *contU); //SobreTodaLaEntrada
uint8_t ConvertBase(int repeat, char rd_ltr);	//Corresponden a la OPER 0=A, 1=C, 2=G, 3=T, 4=N, 5=NN, 6=NNN, 7=sNNNN
char DesconvertBase(uint8_t Ba, int repeat);
char Index2Base(uint8_t Index);

//****************************MAP
void DeConstructOptimizedMap(uint8_t *MAPA , uint64_t  *MAPAO, int lRef, int LenOpMap);
void DeConstructOptimizedMap8(uint8_t *MAPA , uint8_t  *MAPAO, int lRef, int LenOpMap);

//**********************************Coding (Compression) ( Inst --> binary coding )

uint8_t TrdBitInst(int counter, uint8_t  rest, uint8_t  *Oper, uint8_t  *BaseRead, uint8_t BaseRef,
		uint16_t *offset, uint16_t lendesc , char strand, int *aux_i);
uint8_t BitsOperR(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, uint16_t lendesc , int *ii);
uint8_t BitsOperF(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, int *ii );
uint8_t BitsBase(uint8_t BRead, uint8_t BRef);
uint8_t BaseIndex(char Base);
uint8_t Offset( uint16_t offset, uint8_t *rest);
uint8_t Preambulo(uint8_t moreFrags, char strand, uint16_t lendesc , uint16_t EdDis);
uint8_t CreateMask8B(uint8_t pos, uint8_t Long);

void Inst2Bin(uint8_t *BinInst,  uint32_t *posBInst, char strand, uint8_t MoreFrags, uint16_t lendesc, uint16_t *Offsets,
		uint8_t *Oper, uint8_t *BaseRead, uint8_t *BaseRef, FILE *outB, FILE *outb2c, FILE *bin, long i, uint16_t EdDis);

//**********************************Decoding - Decompression (Binary --> Inst)
int CountDelBases(uint16_t Oper);
char DNAComplement(char Base);
int CountInsBases(uint16_t Oper);
char Myitoa(uint8_t val);
void Bin2Inst(uint8_t *BinInst,  uint32_t *posBInst,  char *strand, uint8_t *MoreFrags, uint16_t *lendesc, uint16_t *_Offsets,
		uint8_t *MoreErr, uint8_t *_Oper, uint8_t *BaseRead, uint8_t *BaseRef, uint32_t Nread, FILE *FInst,
		uint8_t *DistReads, uint16_t *CountDels, uint16_t *CountIns);
void Bin2Preambulo(uint8_t FstByte, uint8_t *moreFrags, char *strand);
uint16_t Bin2Offset(uint8_t *BinInst, uint32_t PosInst);
void Byte3Bin2Inst(uint8_t Byte3,uint8_t *MoreErr, uint8_t *oper, uint8_t *DistReads);
uint64_t GetNextPos(uint64_t PrevPos, uint8_t *Map, long lenRef);

//**********************************Builder
void ComplementRead(char *Read, int length);
void ReverseRead(char *Read, int length);
int BuildRead(char strand, uint16_t NDel, uint16_t NIns,  uint16_t LenRead, uint32_t MapPos, char *ref1, long NRead, uint16_t *Of,
		uint16_t lendesc,  uint8_t *Op, uint8_t *Distance,	char *Read /*Olg olg,	uint8_t NInst Instr *MInst, , int lenR, uint8_t Match,  uint8_t *Ba*/);
void CalculateRefsF(char *RefB, uint32_t MapPos, char *ref1, uint16_t *Of, uint16_t lendesc,  uint8_t *Op);


//**********************************Quicksorting
void optimizedQuickSort( uint32_t *A, int low, int high);
void optimizedQuickSortIndex( uint32_t *Pos, int low, int high, uint32_t *Indexes);
int randomizedPartitionInd ( uint32_t *arr, int low, int high, uint32_t *Indexes);
int randomizedPartition ( uint32_t *arr, int low, int high);
void insertionSort( uint32_t *arr, int low, int high);
void insertionSortInd( uint32_t *arr, int low, int high, uint32_t *Indexes);
void optimizeMap(uint8_t *MAPA , uint64_t  *MAPAO, int lRef, int *LenOpMap);
void optimizeMap8(uint8_t *MAPA , uint8_t  *MAPAO, int lRef, int *LenOpMap);


//**********************************File manipulation
void PrintBinary2(unsigned int num, int nbits);
void PrintByte2File(FILE *outB, uint8_t num, int nbits);


int main (){
	
	srand(time(NULL));
	clock_t begin, end;
	begin = clock();
	remove("PreambFile.txt");remove("WeirdFile.txt"); remove("Verifica.txt"); remove ("DecompressedBuilder.txt");

	char Reads[NUM_READSM+3][500];//[READ_LENGTH+1] ; 	//Solo por pruebas del builder ... MUST REORDER THEM !
	int longRead[NUM_READSM+3];			//Solo por pruebas no existe en la version final

/**************************************************CARACTERIZADOR***********************************************************/
/**************************************************END OF CARACTERIZADOR***********************************************************/

/***************************************************CLASIFICADOR**********************************************************/
/***************************************************END OF CLASIFICADOR**********************************************************/

/**************************************************ALINEADOR ******************************************************************/
	//Aligner should point to an specific reference or several chormosomes

/**************************************************END OF ALIGNER ******************************************************************/

	//**************************************REFERENCE READER***********************************************/
		char *ref, *PrevRead , *RefLocation;
		int lenRef; 	//ReferenceLength
		RefLocation="./SamData/lambda_virus.fa";
		ref= load_reference(RefLocation, &lenRef); //LenRef es la length + 1 position used for the zero
		//printf("longitud de la ref %i deberia ser 48502 y es %i \n", strlen(ref), lenRef-1);//En SAM la ref empieza en 1, por lo cual lenREF tiene 1 De MAS


		//***************************************DYNAMIC ARRAYS**************************************************************/
		/**********************************************Dynamic Arrays for the Instructions**************************************************/
			//************PRELUDE
			uint32_t *Position=(uint32_t *)malloc(NUM_READSM*sizeof(uint32_t));    //Mapping position of each read in the reference
			if (Position==NULL) printf ("Not enough memory for Position");

			uint8_t *Match=(uint8_t *)malloc(NUM_READSM*sizeof(uint8_t));      		//Match==Strand, kind of transformation applied before mapping
			if (Match==NULL) printf ("Not enough memory for MAtch");

			uint16_t *CaEr=(uint16_t *)malloc(NUM_READSM*sizeof(uint16_t));    		//Quantity of errors in each read
			if (CaEr==NULL) printf ("Not enough memory for CaEr");

			uint16_t *EdDis=(uint16_t *)malloc(NUM_READSM*sizeof(uint16_t));    	//EditDistance in each read
			if (EdDis==NULL) printf ("Not enough memory for CaEr");

			uint16_t *AdBa=(uint16_t *)malloc(NUM_READSM*sizeof(uint16_t));        //???Cantdad de Bases Adicionales de Este Read, para preveer reads de tam variable
			if (AdBa==NULL) printf ("Not enough memory for AdBa");				   //USO NO IMPLEMENTADO
			//Por ahora en un experimento inicial se le asignara la longitud del READ, porque es calculable al alinear.


			typedef  uint16_t* ptr16; //2D Structure for the offsets, Reads are in rows, Offsets in columns
			ptr16 *Offsets;
			Offsets=(ptr16 *)malloc(NUM_READSM*sizeof(ptr16));
			if (Offsets==NULL) printf ("Not enough memory for Offsets");

			//************OP
			typedef  uint8_t* ptr8; //Por FILAS los reads, por columnas las bases que van en eel read
			ptr8  *OPsBase;
			OPsBase=(ptr8 *)malloc(NUM_READSM*sizeof(ptr8));  //La base de la operacion. Arreglo de punteros
			if (OPsBase==NULL) printf ("Not enough memory for OPsBase");

			//************REF
			typedef  uint8_t* ptr8; //Por FILAS los reads, por columnas los OPs
			ptr8  *BaseRef;
			BaseRef=(ptr8 *)malloc(NUM_READSM*sizeof(ptr8));  //La base de la REF. Arreglo de punteros
			if (BaseRef==NULL) printf ("Not enough memory for BaseRef");

			typedef  unsigned char* ptrC;
			ptrC *OPsCode;
			OPsCode=(ptrC *)malloc(NUM_READSM*sizeof(ptrC));  //La LEtra codigo de la operacion. Arreglo depunteros
			if (OPsCode==NULL) printf ("Not enough memory for OPsCode");
			//Pudiera ser mas eficiente , se desperdician muchos bits

			typedef char* ptrCh;
			ptrCh *UnMapped;
			UnMapped=(ptrCh*)malloc(NUM_READSM*sizeof(ptrCh));  //La LEtra codigo de la operacion. Arreglo depunteros
			if (UnMapped==NULL) printf ("Not enough memory for UnMapped");

			uint8_t *MAPA;
			if (MAPA=(uint8_t *)calloc(lenRef+1,sizeof(uint8_t))) {printf (" Creado mapa ");}else printf (" No se puede crear mapa ");

		int NReads=NUM_READSM+1; //10001;  //GET it from the aligner NUM_READSM,
		int NTErrors=MAX_ERROR_READ+1;//25;    //**** THIS MUST CHANGE*****// Max Number of errors allowed per read


/**************************************************SAM READER********************************************************************/
	//**********************************************Reading the SAM/BAM File **************************************************/
		Olg olg;		//Object for a single read info
		olg.EditDist = 0; 	olg.lendesc=0;
		olg.lonRead=0; 	  		olg.NDel=0;       	olg.NIns=0;

		Instr *MiInstr; //Object for the instructions I get from the olg Object
		int len; 		//Length of the resulting instruction (of each read)
		long int contU=0, contNR=0; //Unmapped reads counter; //contandor de no reconstruidos
		long i=0, PrevTamRead, UC=0;			//Reads Counter (all of them)

		//Get off this const when read it as input
		char const *fname= "./SamData/eg1.sam";     //Needed since in this 1st version data comes from a Sam File instead of from the aligner
		samfile_t *in = samopen(fname, "r", 0);  //This fname should be an input variable
		if (in == 0) {
			fprintf(stderr, "Fail to open SAM file %s\n", fname);
			exit(1);
		}else printf ("SAM File opened \n");
		bam1_t *b=bam_init1();		//bam file variable, optimizing the process

		FILE *out=fopen("Inst_output.txt", "w");    //File with the instruction data after being parsed by this program

		while (samread(in, b) > 0) { 			//only operate on the alignments for this chromosome in the header
			///************************************OUPUT FROM THE SAM FILE (in olg Objects)************************************************///

			olg=MIparse_input(b, ref, lenRef,&contU);  //Parses all the SAM file in the input

			i++; //Reads Counter, Zero postiion remains empty
			olg.NREAD=i;							 //Id of the read (A correlative number)
			fprintf(out, "\n READ %ld ", i);				fprintf(out," LOC %i ", olg.loc);
			fprintf(out,"  descriptor %s", olg.descriptor);
			//fprintf(out,"  Prefijo %s", olg.prefix);	//fprintf(out,"  Sufijo %s", olg.suffix); //Useless
			fprintf(out,"  Read %s",	olg.read);		fprintf(out,"  STRAND %c ", olg.strand);
			fprintf(out,"  #Ins %i, #Dels %i \n", olg.NIns, olg.NDel); 		fflush(out);
			fprintf(out,"  LongREAD %i \n", olg.lonRead); 		fflush(out);

			/***********************Move the DAta to Dynamic Arrays ***********************/
			Position[i]=olg.loc;    //Posicion de MAPEO en REF
			MAPA[olg.loc]=1;		//Building the BitsMap
			Match[i]=olg.strand;

			MiInstr=parse_line(olg, &len, ref); //Descompone el descriptor almacenado en el olg, generando array de instrucciones en las estructuras de memoria dinamica
			olg.lendesc=len;               //Longitud de la instrucion segun el descriptor.
			CaEr[i]=olg.lendesc;     //Cantdad de Errores (EditDist) de Este Read,
			EdDis[i]=olg.EditDist;
			AdBa[i]=olg.lonRead;       //Debe definirse claramente mas tarde, Por el momento corresponde a la cantidad de bases por si se manejan read de L variable

			fprintf(out," Longitud calculada de la Instruccion MD %i , distEdicion %i \n ", olg.lendesc, olg.EditDist); fflush(out);

			strcpy(Reads[i],olg.read) ; 	//Solo por pruebas.
			longRead[i]=olg.lonRead; 		//Solo por pruebas.


			//BE SURE THAT CaEr[i] < MAX_ERROR_READ
			//Ask for memory to store each instruction error data
			Offsets[i]=(uint16_t *)malloc(CaEr[i]*sizeof(uint16_t));
			if (Offsets[i]==NULL) printf ("Not enough memory for Offsets %ld \n",i);

			OPsCode[i]=(unsigned char *)malloc(CaEr[i]*sizeof(unsigned char)); //Codigo de la operacion aplicable a este error
			if (OPsCode[i]==NULL) printf ("Not enough memory for OpsCode %ld \n",i);

			OPsBase[i]=(uint8_t *)malloc(CaEr[i]*sizeof(uint8_t));
			if (OPsBase[i]==NULL) printf ("Not enough memory for OpsBase %ld \n",i);

			BaseRef[i]=(uint8_t *)malloc(CaEr[i]*sizeof(uint8_t));
			if (BaseRef[i]==NULL) printf ("Not enough memory for BaseRef %ld \n",i);

			for (int w=0;w<olg.lendesc; w++){   //Each error in the descriptor is saved and sored
				if (DEBUG){
					fprintf(out,"  Repite %i ",	MiInstr[w].repeat); 			fprintf(out,"  Offset %i ", MiInstr[w].offset);
					fprintf(out,"  Oper %c ", MiInstr[w].type);  				fprintf(out,"  LetRef %c ", MiInstr[w].ref_ltr);
					fprintf(out,"  LetRead %c \n", MiInstr[w].rd_ltr);  		fflush(out);
				}


				//Put Instructions in Dynamic Arrays
				Offsets[i][w]=MiInstr[w].offset;
				OPsCode[i][w]=MiInstr[w].type;
				OPsBase[i][w]=ConvertBase(MiInstr[w].repeat,MiInstr[w].rd_ltr);	//0=A, 1=C, 2=G, 3=T, 4=N, otherwise
				BaseRef[i][w]=ConvertBase(1,MiInstr[w].ref_ltr);

			}//endfor_w=0

			//Initial Checking oh the read reconstrution
			if ((olg.loc!=0)&&(!VerificRead(olg,MiInstr,ref,lenRef, Match[i],Offsets[i],OPsCode[i],OPsBase[i])==0)){
				contNR++; //Wrong Reads Counter
			}

	    }//WHILE (samread(in, b)

		if (DEBUG) fprintf (out,"TOTAL READS %ld de los cuales Mal Reconstruidos %ld NO MAPEADOS %ld \n", i,contNR, contU ); fflush(stdout);

	    if (out) fclose(out);
	    bam_destroy1(b);
	    samclose(in);
	    if (ref!=NULL) free(ref);

/**************************************************END OF SAM READER********************************************************************/

/***********************************************CREATING AND OPTIMIZING THE REFERENCE MAP ***********************************************/
		/*Generacion del mapa optimizado*/ //This is, the bit MAP parallel to the reference//
	    MAPA[0]=0;
		uint64_t *MAPAO;
		//BE SURE THAT lenRef < REF_LENGTH
		if (MAPAO=(uint64_t *)calloc((int)((lenRef+1)/64),sizeof(uint64_t))) {printf ("Creado mapa Opt_64\n");}else printf ("No se puede crear mapa Opt_64\n");
		int LenOpMap;  //All reference LEngth must BE LONG!!!, this and lenREf

		optimizeMap(MAPA , MAPAO, lenRef, &LenOpMap);
		printf ("\n Longitud del mapa optimizado de u64 %i bytes y Ref %i ", LenOpMap+1, lenRef);


		uint8_t *MAPAO_8;
		if (MAPAO_8=(uint8_t *)calloc((int)((lenRef+1)/8),sizeof(uint8_t))) {printf ("Creado mapa Opt_8\n");}else printf ("No se puede crear mapa Opt_8\n");
		int LenOpMap8;  //All reference LEngth must BE LONG!!!, this and lenREf

		optimizeMap8(MAPA , MAPAO_8, lenRef, &LenOpMap8);
		printf ("\n Longitud del mapa optimizado de u7 %i bytes y Ref %i ", LenOpMap+1, lenRef);

		FILE *map; 	/*Escritura del mapa en archivos*/
		if (DEBUG){
			if (map=fopen("MAPA.txt","w")){
			   	for  (long y=1;y<=lenRef;y++) fprintf (map,"%i", MAPA[y]);
			   	fclose(map);
			}
					//if free

			if (map=fopen("MAPA_Opt.txt","w")){
			   	for  (long y=1;y<=lenRef;y++) fprintf (map,"%lu\n", MAPAO[y]);
			   	fclose(map);
				}
		}


		//if (MAPA) free (MAPA);
		//if (MAPAO) free (MAPAO);
		//if (MAPARec) free (MAPARec);
/*********************************************END OF CREATING AND OPTIMIZING THE REFERENCE MAP *********************************************/

/*************************************************SORTER*********************************************************************/

	//THIS OUTPUT SHOULD SEPARATE MAPPED AND UNMAPPED READS IN DIFFERENT ARRAYS
	long tope = i;

	FILE *srt; 	/*Write the original position MAP to File*/  //No order at all, just original reads position
	//Remember position 0 is non valid
    if(DEBUG){
    	if (srt=fopen("UnsortedPositions.txt","w")){
    	  	for  (int y=1;y<=tope;y++) fprintf (srt,"%i\n", Position[y]);
    	   	fclose(srt);
    	}
    }


	//Reserve and Initialize memory for the SORTED index map //In this map are all the contiguos MApPositions
	uint32_t *IndexPos=(uint32_t *)malloc(NUM_READSM*sizeof(uint32_t));    //Mapping position of each read in the reference
	if (IndexPos==NULL) printf ("Not enough memory for IndexPos");
	   for (long g=0; g<=tope; g++) IndexPos[g]=g;

	//Creo un auxiliar solo mara no dañar actualmente el valor de position, luego se deberia eliminar
	uint32_t *AuxPos=(uint32_t *)malloc(NUM_READSM*sizeof(uint32_t));    //Mapping position of each read in the reference
	Position[0]=0;
	if (AuxPos==NULL) printf ("Not enough memory for AuxPos");
	for (long g=0; g<=tope; g++) AuxPos[g]=Position[g];

	   //Call Quicksort , //Phase two, at the time separate mapped from unmapped
	optimizedQuickSortIndex(AuxPos, 1, tope, IndexPos); //optimizedQuickSort(Position, 0, i-1);

	//Verifico que no hayan indices repetidos Could be eliminated in final version
	if (DEBUG) {
		uint32_t *Rep=(uint32_t *)malloc(NUM_READSM*sizeof(uint32_t));
		for (long g=1; g<=tope; g++) Rep[g]=0;
		for (long g=1; g<=tope; g++) {
			if (Rep[IndexPos[g]]==1) {
				printf ("****indice Repetido****");
			}else Rep[IndexPos[g]]=1;
		}
		if (Rep) free(Rep);

		if (srt=fopen("SortedPositions.txt","w")){
		   	for  (int y=1;y<=tope;y++) fprintf (srt,"%i\n", AuxPos[y]);
		   	fclose(srt);
		}
		//AuxPos se puede liberar ya que solo se usa para verificar correcto ordenamiento

		if (srt=fopen("SORTEDIndexes.txt","w")){
		   	for  (int y=1;y<=tope;y++) fprintf (srt,"%i\n", IndexPos[y]);
		   	fclose(srt);
		}

	}

	//Separar mapped de unmmapped, and calculate star of mapped reads
	uint8_t aux=1;
	uint32_t startMapped=0,y=1;
	uint32_t *UnMReIndex=(uint32_t *)malloc(NUM_READSM*sizeof(uint32_t)); //Si se tiene desde el alineador se puede pedir solo la memoria para la cant de unmmaped reads
	while ((aux==1)&&(y<=tope)) {

		if (Position[IndexPos[y]]==0){		//armar unmapped array.
			UnMReIndex[y]=IndexPos[y];
		}else{
			aux=0;
			startMapped=y;//Detectar inicio de MAPPED READS.
		}
		y++;
	}
	printf ("\nLos MAPPED READS INICIAN EN %i: ", startMapped);

	if (DEBUG) {
		FILE *UR; 	/*Escritura en archivos  de reads no mappeado , for developping purposes only*/ //In this map are all the contiguos MApPositions
		if (UR=fopen("UnmappedReads.txt","w")){
			for  (int y=1;y<startMapped;y++) fprintf (UR,"%s\n", Reads[UnMReIndex[y]]);
			fclose(UR);
		}
	}
    //Enviar al ORCOM  los NO MAPEADOS Y LUEGO LIBERAR
	/*for(int z=0; z<NUM_READSM;z++){ //EY SOLO HASTA EL CONTADOR DE READS NO MAPEADOS
		if (UnMapped[z]) free (UnMapped[z]);
	};*/

	/*ENTONCS EL ARREGLO DE INDEX POS SE PUEDE PICAR EN DOS PARTES desde [1..startMapped) los Reads no mapeados [startMapped, tope] los mapeadso
	 * Ahora , yo para optimizar el computo podria pedir que este ordenamiento venga ya del alineador listo y trasladar este computo par aalla.
	 *
	//ENTONES SE PUEDE APROVECHAR ESTO PARA QUE BUILDER Y RESTO TRABAJEN EXLUSIVAMENTE SOBRE LOS MAPPED


	/* Liberar aca y luego despues del BUILDER EN LAS EST CREADAS PARA TAL FIN
	//**Liberando la memoria usada para estructurar las instrucciones, NO YET e lbuilder depende esto
	if (Position) free(Position);
	*/


	//Distribute output to the different cores and UNMAPPED TO THE GPU

/*************************************************END OF SORTER*********************************************************************/

/***************************************************COMPRESSOR****************************************************************/
	uint8_t *BinInst; //In the BinInst array we will store all of the 3 bytes (or less) to describe each read
	if (BinInst=(uint8_t *)calloc((NReads*NTErrors*BYTES_PER_ERROR),sizeof(uint8_t))) {printf ("\n");}else printf ("No se puede crear BinInst\n");
	//0.Reservar la Memoria para los 3 enteros: Cuantos:[Cantidad de Reads]+[CantidadTotalErroresx2] tipo: unint8
	//Estrategia: CantReads COnocida. , porque viene del alineamiento, ojo con cant de operaciones usadas.

	uint32_t posBInst=0;//Current position in the BinInst array, will go to 1 soon
	uint32_t ReadLength;//Comes from the aligner

	//Instruction Builder CODER

	FILE *out2=fopen("ReadsOnly.txt", "w");	    //File with the Reads Only as extracted from the Fastq File
	FILE *outb2c=fopen("Byte2UInt.txt", "w");   //File with the calculated instruction in Unsigned Int Representation
	FILE *outbin=fopen("Byte2Binary.bin", "wb");//Binary File with the calculated instruction in byte
	FILE *outB=fopen("Bytes.txt", "w"); 		//Bytes Stream File
	uint8_t MoreFrags=0;
	long AuxInd;

	for (long i=startMapped ; i<=tope; i++ ){  //Only for mapped reads, //Position[]==0 is the same than unmapped
		AuxInd=IndexPos[i];
		fprintf(out2,"%s\n",Reads[AuxInd]);//Prints the only read file //EYE ON IT

		if ((i<tope)&&(Position[AuxInd]!=0)&&(Position[AuxInd]==Position[IndexPos[i+1]])) MoreFrags=1; //Always check the NEXT mapping postition

		//*************CODING to bin the current read
		if ((Match[AuxInd]!='U')) //*************CODING to bin the current read, //Only for mapped reads
			Inst2Bin(BinInst, &posBInst, Match[AuxInd], MoreFrags, CaEr[AuxInd], Offsets[AuxInd],OPsCode[AuxInd], OPsBase[AuxInd], BaseRef[AuxInd], outB, outb2c, outbin, AuxInd, EdDis[AuxInd]);

		MoreFrags=0;
    }// // for long startM
	if (DEBUG){
		printf ("A Main.- The BinInst [1-4]  is %u %u %u %u \n",BinInst[1],BinInst[2],BinInst[3],BinInst[4]);
		printf ("A Main.- The BinInst [10-100-200]  is %u %u %u |\n",BinInst[10],BinInst[100],BinInst[200]);
		printf ("A Main.- The BinInst [10k-25k-53k]  is %u %u %u |\n",BinInst[10000],BinInst[25000],BinInst[53000]);

	}

	if (outB) fclose(outB);
	if (out2) fclose(out2);
	if (outb2c) fclose(outb2c);
	if (outbin) fclose(outbin);

	//if (Match)    free (Match); //Uncomment in the final fphase since decompressor depends on this, Same as Position
	if (CaEr)     free (CaEr);
	if (EdDis)    free (EdDis);
	if (AdBa)	  free (AdBa);

	//free (Moerr, mediante un ciclo);
	for(int z=0; z<NUM_READSM;z++){ //HASTA LA CANTIDAD DE READS MAPEADOS OJO
		if (UnMapped[z]) free (UnMapped[z]);
		if (Offsets[z]) free(Offsets[z]);
		if (OPsCode[z]) free(OPsCode[z]);
		if (OPsBase[z]) free(OPsBase[z]);
	};

	/***************************************************END OF HIGH LEVEL COMPRESSOR****************************************************************/


/************************************************LOW LEVEL COMPRESSOR (pbzip2 HLI)********************************************************************/
	//PROBLEMAS : NO CORRE EN PARALELO. GENERA DIRECTAMENTE EL ARCHIVO COMPRIMIDO (buscar que genere bloque de memoria), concatenar, empaquetar y vaciar archivo.
	//LA 2da corrida el BInInst tiene valores diferentes: OJO.

	FILE *fpout=fopen("LowLevelData.bz2", "wb");//Binary File with the calculated instruction in byte;
	/*Comprimir en Bzip : tamaño de los reads (fixed) , mapa, instrucciones, ... Unmmaped reads Orcom compressed ?*/
	int bzerror = BZ_OK;
	BZFILE *bfp = BZ2_bzWriteOpen(&bzerror, fpout, 9, 0, 30);

	/* BZFILE *BZ2_bzWriteOpen( int *bzerror,   FILE *f,     int blockSize100k,  int verbosity,            int workFactor );

- blockSize100k specifies the block size to be used for compression. This number x 100.000 bytes
- verbosity should be set to a number between 0 and 4 inclusive.   0 is silent, and greater numbers give increasingly verbose monitoring/debugging output
- workFactor controls how the compression phase behaves when presented with worst case,   highly repetitive, input data.
  If compression runs into difficulties caused by repetitive data,the library switches from the standard sorting algorithm to a fallback algorithm.
  The fallback is slower than the standard algorithm by perhaps a factor of three, but always behaves reasonably, no matter how bad the input
	 * */

	if (bzerror != BZ_OK)
	{
	    BZ2_bzWriteClose(&bzerror, bfp, 0, NULL, NULL);
	    //fclose(fpin); si tuviesemos un archivo de entrada que cerrar
	    fclose(fpout);
	    return 1;
	}



	//1 Read RefLocation= --> 32 bytes. (location). Dificultad: Manejar tamaño
	// Finish it
	int leng = 32;
	char buf32[32]={0};//=itoa(lengthReads,10);
	if (strlen(RefLocation)<32) strcpy(buf32, RefLocation);
	//snprintf(buf, leng+1, "%d", lengthReads);

    	BZ2_bzWrite(&bzerror, bfp, buf32, leng);  //Variable de error, en el archivo que va, la data del buffer, longitud.
    	if (bzerror == BZ_IO_ERROR){
        	printf("bz-io-error detected\n");
    	}


	//2.-Read LengthRef = 64 bits--> 8 bytes.. up to 99.999.999
	uint32_t lengthRef=lenRef; //8 bytes
	char buf8[8]={0};//=itoa(lengthReads,10);
	/*int */leng=8;
	snprintf(buf8, leng+1, "%d", lengthRef);
	BZ2_bzWrite(&bzerror, bfp, buf8, leng);  //Variable de error, en el archivo que va, la data del buffer, longitud.
	if (bzerror == BZ_IO_ERROR){
    	printf("bz-io-error detected\n");
	}

/*
	//3.-Length MapOP8 --> 8 bytes.. up to 99.999.999
	uint32_t lengthMap8O=68439; //8 bytes
		 buf8[8]={0};//=itoa(lengthReads,10);

		snprintf(buf8, leng+1, "%d", lengthMap8O);
		BZ2_bzWrite(&bzerror, bfp, buf8, leng);  //Variable de error, en el archivo que va, la data del buffer, longitud.
		if (bzerror == BZ_IO_ERROR){
	    	printf("bz-io-error detected\n");
		}
*/
	//3.-LEngth Bin Inst --> 8 bytes.. up to 99.999.999
	uint32_t lengthBInst=posBInst; //8 bytes
	buf8[8]={0};//=itoa(lengthReads,10);
	snprintf(buf8, leng+1, "%d", lengthBInst);
	BZ2_bzWrite(&bzerror, bfp, buf8, leng);  //Variable de error, en el archivo que va, la data del buffer, longitud.
	if (bzerror == BZ_IO_ERROR){
	   	printf("bz-io-error detected\n");
	}

	//4.- GUARDAR LONG DE LOS READS (EXperiment with any 32 integer), convertirlo a atoi
	uint32_t lengthReads=8439; //32 bits //Como separarlo en 4 cadenas de 8 bits?
	/*int leng=1;
	if ((lengthReads>=10)&&(lengthReads<100)) leng=2;
		else{
			if ((lengthReads<1000)) leng=3;
				else{
					if (lengthReads<10000) leng=4;
						else leng =5;
				}
			}*/
	/*int*/ leng = 5;
	char buf[5]={0};//=itoa(lengthReads,10);
	snprintf(buf, leng+1, "%d", lengthReads);

    	BZ2_bzWrite(&bzerror, bfp, buf, leng);  //Variable de error, en el archivo que va, la data del buffer, longitud.
    	if (bzerror == BZ_IO_ERROR){
        	printf("bz-io-error detected\n");
    	}

	//5 Save Map (8 bits version OJO)
	// En 64 bits caben 8 caracteres simples. En 100 k caben hasta 1024 bytes o sea 128 enteros  de 64 bits OJO
    char *PtrInst; int BufCapacity=899000;//1024;//Should be 900.000 bytes
    long IndexInst=0;

	PtrInst=(char *)MAPAO_8;

	if (DEBUG) {
    	printf ("A.- W The PtrInst [1-3]  is %u %u %u \n",PtrInst[10],PtrInst[100],PtrInst[200]);
    	printf ("A.- W The MAPAO_8 [1-3]  is %u %u %u |\n",MAPAO_8[10],MAPAO_8[100],MAPAO_8[200]);

    	printf ( "IndexInst %i Diferencia %i \n" ,IndexInst, (int)(LenOpMap8+1)-BufCapacity );
	}
	while ((bzerror != BZ_IO_ERROR)&&(IndexInst<=((int)(LenOpMap8+1)-BufCapacity))){
        printf ("ICI\n");
		if ((DEBUG) && (IndexInst<BufCapacity)){ //Only HERE for Dubugging
			printf ("B.- W The PtrInst [1-3]  is %u %u %u \n",PtrInst[2000],PtrInst[5000],PtrInst[5500]);
			printf ("B.- W The MAPAO_8 [1-3]  is %u %u %u |\n",MAPAO_8[2000],MAPAO_8[5000],MAPAO_8[5500]);
		}

	    BZ2_bzWrite(&bzerror, bfp, PtrInst, BufCapacity);  //Variable de error, en el archivo que va, la data del buffer, longitud.

	    PtrInst=PtrInst+BufCapacity;
	    IndexInst=IndexInst+BufCapacity;

	    if (bzerror == BZ_IO_ERROR){        		printf("bz-io-error detected\n");    		}
	}//end while

	long diff = (int)(LenOpMap8+1)-IndexInst; 	    //procesar el resto de datos en el bbuffer

	if (diff>0) BZ2_bzWrite(&bzerror, bfp, PtrInst, diff);  //Variable de error, en el archivo que va, la data del buffer, longitud.

	if (bzerror == BZ_IO_ERROR){	     			printf("bz-io-error detected\n");	    	}

	//6.--- Save instructions, 	    	//Mis datos son arreglo de 8 bts, lo paso directo a byte

    IndexInst=0;
	PtrInst=(char *)BinInst;
	if (DEBUG){
		printf ("A.- The BinInst [1-4]  is %u %u %u %u \n",BinInst[1],BinInst[2],BinInst[3],BinInst[4]);
		printf ("A.- W The PtrInst [10-100-200]  is %u %u %u \n",PtrInst[10],PtrInst[100],PtrInst[200]);
		printf ("A.- W The BinInst [10-100-200]  is %u %u %u |\n",BinInst[10],BinInst[100],BinInst[200]);
	}

	while ((bzerror != BZ_IO_ERROR)&&(IndexInst<= (int)(posBInst+1)-BufCapacity)){

		if ((DEBUG)&&(IndexInst<1024)){ //Only HERE for Dubugging
			printf ("B.- W The PtrInst [10k-25k-53k]  is %u %u %u \n",PtrInst[10000],PtrInst[25000],PtrInst[53000]);
			printf ("B.- W The BinInst [10k-25k-53k]  is %u %u %u |\n",BinInst[10000],BinInst[25000],BinInst[53000]);
		}

	    BZ2_bzWrite(&bzerror, bfp, PtrInst, BufCapacity);  //Variable de error, en el archivo que va, la data del buffer, longitud.

	    PtrInst=PtrInst+BufCapacity;
	    IndexInst=IndexInst+BufCapacity;


	    if (bzerror == BZ_IO_ERROR){        		printf("bz-io-error detected\n");    break;		}
	    //memset(PtrInst, 0, BufCapacity); //Re assign into final variable
	}//end while
	diff = (posBInst+1)-IndexInst; 	    //procesar el resto de datos en el bbuffer

	if (diff>0) BZ2_bzWrite(&bzerror, bfp, PtrInst, diff);  //Variable de error, en el archivo que va, la data del buffer, longitud.

	if (bzerror == BZ_IO_ERROR){	     			printf("bz-io-error detected\n");	  		 	}

	//7.-Release all
 	BZ2_bzWriteClose(&bzerror, bfp, 0, NULL, NULL);
 	//free poninters
	fclose(fpout);



	/************************************************END OF LOW LEVEL COMPRESSOR (pbzip2 HLI)********************************************************************/


	/************************************************LOW LEVEL COMPRESSOR (pbzip2 LLI)********************************************************************/

	//NOTAS DE DESARROLLO: En lo consiguiente, se parte de un ejemplo por lo cual, replicare toda la linea 2 veces, la primera aparicion, comentada generalemnte
    //correspondera a la sintaxis original segun aparecia en el ejemplo. todas las corerspondientes debajo de esta seran intentos o versiones mias
    //SI una linea aparece solo 1 vez coincide en ser igual en mi version qu een la original, sin cambios

	//Ejemplos:
	    // https://stackoverflow.com/questions/9577735/how-compress-data-in-memory-buffer-by-using-libbz2-library-in-c-program
	   //https://stackoverflow.com/questions/13065023/using-bzip2-low-level-routines-to-compress-chunks-of-data

	/*Totalidad de Elementos a comprimir:
	 	 * leng (direccion de la ref): 32 bytes
		 * +lengthRef (long de la ref): 8 bytes
		 * lengthBInst (PosbInst ) (long del arreglo a comprimir): (longitud del arreglo de instrucciones): 8 bytes.
		 * lengthReads (longitud de los reads): 4 bytes
		 * 		NO NECESITA GUARDARSE, es la misma de las referenciaLenOpMap8 (longitud del MAPA):
		 * + PosbInst* 1byte each
		 *
		 * Entonces totalidad de bytes a comprimir: 32 bytes+8bytes+8bytes+4 bytes+ (lengthBInst) bytes<tantos bytes como elementos en el areglo de instrucciones>)
		 * Entones: 32 bytes+8bytes+8bytes+4 bytes+ (lengthBInst) = (52 + lengthBInst)Bytes
		 * MAximo de bytes a comprimir en UNA SOLA ITERACION (bzip defined):
		 * */

	//Duda : Donde , como y sobre quien se hace la reserva de memoria? una variable x creada. que taaño ?
	//Usar un malloc para reservar memoria para el maximo de lementos a comprimir en una iteracion
	//En bytes: cuantos elementos se han de comprimir: = + + + =


//Definir el valor de work factor
//https://stackoverflow.com/questions/13065023/using-bzip2-low-level-routines-to-compress-chunks-of-data

uint64_t TotalData=0UL ; //Total amount of data to compress in bytes (# of bytes required)
bz_stream bzStream;
int64_t bzBytesWritten = 0UL;
uint64_t cumulativeBytesWritten = 0ULL;
size_t myBufferLength = 0;
uint64_t TotalBytestoCompress = (52 + lengthBInst);//ajgs /***OJO FALTA SUMARLE LA CANTIDAD DE BYTES DEL MAPA****///

unsigned char bzBuffer[BZIP2_BUFFER_MAX_LENGTH] = {0}; //darle valor previo a esta const
unsigned char myBuffer[UNCOMPRESSED_MAX_LENGTH] = {0}; //darle valor previo a esta const
FILE *outLLIComp=fopen("LLICOmpression.som", "b");


//.- Prepares for compression:  allocate a bz_stream
	// int BZ2_bzCompressInit ( bz_stream *strm, int blockSize100k, int verbosity, int workFactor );
	//dedicarle un tiempo a decidir el valor de work facto como explica la pagina 13 del manual

	 // initialize bzStream
	 // By Me August 2017
	  bzStream.next_in = NULL;
	  bzStream.avail_in = 0U; // [atencion en que valor iniciar este ? ]
	  bzStream.avail_out = 0U; // [atencion en que valor iniciar este ? ]
   	  bzStream.bzalloc = NULL;
	  bzStream.bzfree = NULL;
	  bzStream.opaque = NULL;
	  int bzError = BZ2_bzCompressInit(&bzStream, 9, 0, 30);

	  if (bzError!=BZ_OK){ //bzError checking...
		  printf ("Severe Error in BzCompressInit, program aborted.");
		  exit(0);
	  }

//1. First Block to compress LLI: RefLocation= --> 32 bytes. (location). Dificultad: Manejar tamaño

	  // compress bytes in myBuffer
	  //bzStream.next_in = myBuffer;
	  /*LET THE DATA BE READY IN buf32 as in the code above*/
	  bzStream.next_in = buf32;//(char *)malloc((TotalBytestoCompress)*sizeof(char)); //ajgs, this should be pointin DTA TO BE COMPRESSED
	  //bzStream.avail_in = myBufferLength;
	  bzStream.avail_in = leng;
	  //bzStream.next_out = bzBuffer;
	  char *bz_Buffer=(char *)malloc((TotalBytestoCompress)*sizeof( char)); //ajgs /*Ojo con total bytes to compress and BZIP2_BUFFER_MAX_LENGTH*/
	  bzStream.next_out= bz_Buffer; //ajgs
	  //bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;
	  bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;//leng; // ajgs

	  //Ahora si comprimo
	  bzError = BZ2_bzCompress(&bzStream, BZ_RUN);

	  if (bzError!=BZ_OK){ //bzError checking...
	  	  printf ("Severe Error in BZ2_bzCompress for RefLocation , program aborted.");
	  	  exit(0);
	  }

	  //bzBytesWritten = sizeof(bzBuffer) - bzStream.avail_out;
	  bzBytesWritten = sizeof(bz_Buffer) - bzStream.avail_out; //may to control the amount of data compressed
	  //cumulativeBytesWritten += bzBytesWritten;
	  cumulativeBytesWritten += bzBytesWritten;



//2. Second Block to compress LLI: LengthRef, Length of the reference
		//2.-Read LengthRef = 64 bits--> 8 bytes.. up to 99.999.999

		//char buf8[8]={0};
		/*int */leng=8;
		//snprintf(buf8, leng+1, "%d", lengthRef);


	  //Based on what I made for section 1
	 // compress bytes in myBuffer

			  bzStream.next_in = buf8;//(char *)malloc((TotalBytestoCompress)*sizeof(char)); //ajgs, this should be pointin DTA TO BE COMPRESSED
			  bzStream.avail_in = leng;
			  //bzStream.next_out = bzBuffer;
			 // bzStream.next_out= bz_Buffer; //ajgs, no se reasigna porque fue asignado en 1
			  //bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;
			  //bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;//leng; // ajgs //ajgs, no se reasigna porque fue asignado en 1

			  //Ahora si comprimo
			  bzError = BZ2_bzCompress(&bzStream, BZ_RUN);

			  if (bzError!=BZ_OK){ //bzError checking...
			  	  printf ("Severe Error in BZ2_bzCompress for RefLocation , program aborted.");
			  	  exit(0);
			  }

			  //bzBytesWritten = sizeof(bzBuffer) - bzStream.avail_out;
			  //bzBytesWritten = sizeof(bz_Buffer) - bzStream.avail_out; //may to control the amount of data compressed
			  //cumulativeBytesWritten += bzBytesWritten;
			  //cumulativeBytesWritten += bzBytesWritten;


//3. Third Block to compress LLI: Length tf BinInst

			  //3.-LEngth Bin Inst --> 8 bytes.. up to 99.999.999
			  	/*uint32_t*/ lengthBInst=posBInst; //8 bytes
			  	buf8[8]={0};//=itoa(lengthReads,10);
			  	snprintf(buf8, leng+1, "%d", lengthBInst);

			  	 //Based on what I made for section 1
			  		 // compress bytes in myBuffer

			  		bzStream.next_in = buf8;//(char *)malloc((TotalBytestoCompress)*sizeof(char)); //ajgs, this should be pointin DTA TO BE COMPRESSED
			  		 bzStream.avail_in = leng;//forzar leng == 8, forzr que este valor entero ocupe siempre 8 bytes
			  				  //bzStream.next_out = bzBuffer;
			  				  //bzStream.next_out= bz_Buffer; //ajgs
			  				  //bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;
			  				  //bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;//leng; // ajgs

			  				  //Ahora si comprimo
			  				  bzError = BZ2_bzCompress(&bzStream, BZ_RUN);

			  				  if (bzError!=BZ_OK){ //bzError checking...
			  				  	  printf ("Severe Error in BZ2_bzCompress for RefLocation , program aborted.");
			  				  	  exit(0);
			  				  }

			  				  //bzBytesWritten = sizeof(bzBuffer) - bzStream.avail_out;
			  				  //bzBytesWritten = sizeof(bz_Buffer) - bzStream.avail_out; //may to control the amount of data compressed
			  				  //cumulativeBytesWritten += bzBytesWritten;
			  				  //cumulativeBytesWritten += bzBytesWritten;



//4. Fourth Block to compress: Length of reads [sequential]


			   /*uint32_t*/ lengthReads=8439; //32 bits //Como separarlo en 4 cadenas de 8 bits?
			  					leng = 5;
			  					/*char*/ buf[5]={0};//=itoa(lengthReads,10);
			  					snprintf(buf, leng+1, "%d", lengthReads);

			  	bzStream.next_in = buf;//(char *)malloc((TotalBytestoCompress)*sizeof(char)); //ajgs, this should be pointin DTA TO BE COMPRESSED
			  	bzStream.avail_in = leng;//forzar leng == 8, forzr que este valor entero ocupe siempre 8 bytes
			  							  				  //bzStream.next_out = bzBuffer;
			  							  				  //bzStream.next_out= bz_Buffer; //ajgs
			  							  				  //bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;
			  							  				  //bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;//leng; // ajgs

			  							  				  //Ahora si comprimo
			  	bzError = BZ2_bzCompress(&bzStream, BZ_RUN);

			  							  				  if (bzError!=BZ_OK){ //bzError checking...
			  							  				  	  printf ("Severe Error in BZ2_bzCompress for RefLocation , program aborted.");
			  							  				  	  exit(0);
			  							  				  }


			  							  				//bzBytesWritten = sizeof(bzBuffer) - bzStream.avail_out;
			  	/*Following sentence might need to be corrected since cz_buffer has chaged a lot in this time.*/
			  	bzBytesWritten = sizeof(bz_Buffer) - bzStream.avail_out; //may to control the amount of data compressed
			  							  							  //cumulativeBytesWritten += bzBytesWritten;
			  	cumulativeBytesWritten += bzBytesWritten;


			  // write compressed data in bzBuffer to standard output
			  	  //fwrite(bzBuffer, 1, bzBytesWritten, stdout);
			  	  fwrite(bz_Buffer, 1, bzBytesWritten, outLLIComp);
			  	  //size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream) writes data from the array pointed to, by ptr to the given stream.
			  	  //El 1 hay q cambiarlo por el tamaño en bytes de cada dato a ser escrito, quizas aca sigue siendo 1, stdout nombre de archivo de salida
			  	  fflush(outLLIComp);
			        //In my case, I don't want to output that data directly, or mabe I do is the data is sincronized
			  	  //and compressed in order.

/*****LAS SIGUIENTES COMPRESIONES SON DIFERENTES PORQUE SON ITERATIVAS; ;TODOS LOS BUFERS Y CONTADORES SON RESETEADOS OK*/

//5. Fifth Block to compress: Map of size [Iterative] 5 Save Map (8 bits version OJO)
	// En 64 bits caben 8 caracteres simples. En 100 k caben hasta 1024 bytes o sea 128 enteros  de 64 bits OJO
	//Posible oportunidad de optimizar segun sugerencia del profe sebas.


			  	// Codigo Guia en la parte de abajo. A continuacion se comprimen los datos en dos fases. la primera en bloques de tamaño masixmo, la segunda  paa el cnjunto restante al final de la compreison
				uint32_t auxPtr= 0, MyBufferLength = 7000; // OJO es este el tamaño correcto ? o //Old int BufCapacity=899000;
				/*char *PtrInst*/; PtrInst=(char *)MAPAO_8;
//Nota ajgs: Para mi aqui hay un ciclo innecesario, no manejo chunks d tamaño discreto. //Main focus by me: a cycle and the rest, as simple as the HLI version


			  	    // read some bytes into myBuffer... // esto es como si se cargara un bloquesote largo para comprimirlo por partes

			  	    bzStream.next_in = PtrInst+auxPtr;         //bzStream.next_in = myBuffer;   //Input//
			  	    bzStream.avail_in = MyBufferLength;    	  //The same. Pero poner la longitud del buffer de entrada

          	  	  	if (bz_Buffer) free(bz_Buffer);           //Free it before using it again

          	  	  	//OJO ACTUALICE O TENGA CUIDADO DE ESTA DIMENSION, //Puede que  no sea una cantidad total sino una parcial
          	  	  	bz_Buffer=(char *)malloc((TotalBytestoCompress)*sizeof(char)); //ajgs :Ojo con total bytes to compress and BZIP2_BUFFER_MAX_LENGTH*/

     	  	      do //Phase II, aqui se comprime por partes ese gran bloque mientras haya dato (cuidar condicion BZ_Ok)
			  	    {      //ajgs: que hace: hay un bloque que se comprime y se va vaciando de una vez al archivo mientras no haya error.
			  	        bzStream.next_out = (char *)bzBuffer;//hay un bloque que se comprime y se va vaciando de una vez al archivo mientras no haya error.
			  	        bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH; //Cuidar este tamaño ver pag 13
			  	        bzError = BZ2_bzCompress(&bzStream, BZ_RUN);

			  	        if (bzError == BZ_IO_ERROR) { printf("bz-io-error detected\n");; break;} // error checking... OJO IMPLEMENTAR, y posiblemente indicar que finalice el ciclo para q no se actualicen los valoes no ejecute el fwrite

                        //actual valor del auxPtr
			  	        bzBytesWritten = sizeof(bzBuffer) - bzStream.avail_out;

			  	        cumulativeBytesWritten += bzBytesWritten;  //Old en mi caso era long IndexInst=0; //los q ya he comprimdo

			  	        fwrite(bzBuffer, 1, bzBytesWritten, stdout); // write compressed data in bzBuffer to standard output
			  	        fflush(stdout);
			  	    }while (/*(bzError == BZ_OK)&&*/(bzError != BZ_STREAM_END)); // &&(IndexInst<=((int)(LenOpMap8+1)-BufCapacity))
//Pregunta : debo tener control de la dimesion de la estructra de datos ? se actualiza automatico el _IN ?




//6. Sixth and last Block to compress: Instructions of size [Iterative]
/*
//aparentemente en next out hay q asignar un puntero con la cantidad de memoria (reservada y vacia) requerida para el dato a comprimir
//aparentemente en avail_out hay que asignar la cantidad de elementos que se colocaran en next_ out, sera igual a next in ?
//Adicionalmente se debe llevar otra variable auxiliar donde se van copiando los elementos que se han comprimido en el paso actual
 //Todo parece indicar que el valor next_in debe apuntar a la posicion de memoria a comprimir
*/


				/*uint32_t*/ auxPtr= 0, MyBufferLength = 7000; // OJO es este el tamaño correcto ? o //Old int BufCapacity=899000;
				/*char *PtrInst*/; PtrInst=(char *)BinInst;
//Nota ajgs: Para mi aqui hay un ciclo innecesario, no manejo chunks d tamaño discreto. //Main focus by me: a cycle and the rest, as simple as the HLI version


			  	    // read some bytes into myBuffer... // esto es como si se cargara un bloquesote largo para comprimirlo por partes

			  	    bzStream.next_in = PtrInst+auxPtr;         //bzStream.next_in = myBuffer;   //Input//
			  	    bzStream.avail_in = MyBufferLength;    	  //The same. Pero poner la longitud del buffer de entrada

          	  	  	if (bz_Buffer) free(bz_Buffer);           //Free it before using it again

          	  	  	//OJO ACTUALICE O TENGA CUIDADO DE ESTA DIMENSION, //Puede que  no sea una cantidad total sino una parcial
          	  	  	bz_Buffer=(char *)malloc((TotalBytestoCompress)*sizeof(char)); //ajgs :Ojo con total bytes to compress and BZIP2_BUFFER_MAX_LENGTH*/

     	  	      do //Phase II, aqui se comprime por partes ese gran bloque mientras haya dato (cuidar condicion BZ_Ok)
			  	    {      //ajgs: que hace: hay un bloque que se comprime y se va vaciando de una vez al archivo mientras no haya error.
			  	        bzStream.next_out = (char *)bzBuffer;//hay un bloque que se comprime y se va vaciando de una vez al archivo mientras no haya error.
			  	        bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH; //Cuidar este tamaño ver pag 13
			  	        bzError = BZ2_bzCompress(&bzStream, BZ_RUN);

			  	        if (bzError == BZ_IO_ERROR) { printf("bz-io-error detected\n");; break;} // error checking... OJO IMPLEMENTAR, y posiblemente indicar que finalice el ciclo para q no se actualicen los valoes no ejecute el fwrite

                        //actual valor del auxPtr
			  	        bzBytesWritten = sizeof(bzBuffer) - bzStream.avail_out;

			  	        cumulativeBytesWritten += bzBytesWritten;  //Old en mi caso era long IndexInst=0; //los q ya he comprimdo

			  	        fwrite(bzBuffer, 1, bzBytesWritten, stdout); // write compressed data in bzBuffer to standard output
			  	        fflush(stdout);
			  	    }while (/*(bzError == BZ_OK)&&*/(bzError != BZ_STREAM_END)); // &&(IndexInst<=((int)(LenOpMap8+1)-BufCapacity))
//Pregunta : debo tener control de la dimesion de la estructra de datos ? se actualiza automatico el _IN ?






//OLD Guide code
/*do{ //Phase I, blocks of maximum and equal size
    // read some bytes into myBuffer...

    bzStream.next_in = myBuffer;
    bzStream.avail_in = myBufferLength;
    bzStream.next_out = bzBuffer;
    bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;
    do
    {
        bzStream.next_out = bzBuffer;
        bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;
        bzError = BZ2_bzCompress(&bzStream, BZ_RUN);

        // error checking...

        //bzBytesWritten = ((unsigned long) bzStream.total_out_hi32 << 32) + bzStream.total_out_lo32;//This one does not works
        bzBytesWritten = sizeof(bzBuffer) - bzStream.avail_out;

        cumulativeBytesWritten += bzBytesWritten;

        // write compressed data in bzBuffer to standard output
        fwrite(bzBuffer, 1, bzBytesWritten, stdout);
        fflush(stdout);
    }while (bzError == BZ_OK);

} while (/* while there is a non-final myBuffer full of discrete chunks left to compress... *///);


// read in the final batch of bytes into myBuffer (with a total byte size of `myBufferLength`...
//Phase II: compress the remaining elements in at the end of the data array
// compress remaining myBufferLength bytes in myBuffer

/*

bzStream.next_in = myBuffer;
bzStream.avail_in = myBufferLength;
bzStream.next_out = bzBuffer;
bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;
do
{
    bzStream.next_out = bzBuffer;
    bzStream.avail_out = BZIP2_BUFFER_MAX_LENGTH;
    bzError = BZ2_bzCompress(&bzStream, (bzStream.avail_in) ? BZ_RUN : BZ_FINISH);

    // bzError error checking...

    // increment cumulativeBytesWritten by `bz_stream` struct `total_out_*` members
     bzBytesWritten = sizeof(bzBuffer) - bzStream.avail_out;

   //bzBytesWritten = ((unsigned long) bzStream.total_out_hi32 << 32) + bzStream.total_out_lo32;
    cumulativeBytesWritten += bzBytesWritten;

    // write compressed data in bzBuffer to standard output
  fwrite(bzBuffer, 1, bzBytesWritten, stdout);
    fflush(stdout);
}while (bzError != BZ_STREAM_END);

*/




// close stream
bzError = BZ2_bzCompressEnd(&bzStream);

// bzError checking [Complete please ajgs!]





	/************************************************END OF LOW LEVEL COMPRESSOR (pbzip2 LowLevInt)********************************************************************/



/********************************************************PACKAGER & OUTPUT WRITER********************************************************/
	//Here and for one time only you should write the read lentgh (0-65535) deberia ser de 2 a la 16 , to the binary file, obviously converted to bin
	/*     //whatever it is
	 //Do not forget to open and close such file fp!
    fwrite((const void*)&ReadLength,sizeof(int),1,fp);
    */

/********************************************************END PACKAGER & OUTPUT WRITER********************************************************/

/****************************************************UNPACKAGER********************************************************************/
/****************************************************END OF UNPACKAGER********************************************************************/


/************************************************LOW LEVEL DECOMPRESSOR (pbzip2 HLI HIGH LEVEL INTERFACE)********************************************************************/

	FILE *fpin=fopen("LowLevelData.bz2", "rb");//Binary File with the calculated instruction in byte;
	//Reservar memoria para el buffer
	bzerror = BZ_OK;
	bfp = BZ2_bzReadOpen(&bzerror, fpin, 0, 0, NULL, 0);
	/* BZFILE *BZ2_bzReadOpen( int *bzerror, FILE *f,        int verbosity, int small,        void *unused, int nUnused ); */

	if (bzerror != BZ_OK)	{
	    BZ2_bzReadClose(&bzerror, bfp);
	    //fclose(fpin); si tuviesemos un archivo de entrada que cerrar
	    fclose(fpin);
	    return 1;
	}

	//Read NameRef= --> 32 bytes.
	//Read LengthRef = 64 bits--> 8 bytes
	//LenOp8 --> 8 bytes
	//LenBintInst --> 8 bytes
	//opcion guardar estos valores o no guardarlos y comprimir/ por  bloque independientes


	//1.- Read RefLocation= --> 32 bytes. (location). Dificultad: Manejar tamaño

	leng = 32;
	/*char*/ buf32[32]={0};//=itoa(lengthReads,10);
	//if (strlen(RefLocation)<32) strcpy(buf32, RefLocation);
	//snprintf(buf, leng+1, "%d", lengthReads);

	BZ2_bzRead(&bzerror, bfp, (char *)buf32, leng);  //Variable de error, en el archivo que va, la data del buffer, longitud.

    if (bzerror == BZ_IO_ERROR){
       	printf("bz-io-error detected\n");
    }
    printf (" Original Location %s , Dec %s  \n", RefLocation, buf32);

	//2.-Read LengthRef = 64 bits--> 8 bytes.. up to 99.999.999
	/*uint32_t*/ lengthRef=0; //8 bytes
	/*char*/ buf8[8]={0};//=itoa(lengthReads,10);
	leng=8;
	BZ2_bzRead(&bzerror, bfp, (char *)buf8, leng);
	if (bzerror == BZ_IO_ERROR){
    	printf("bz-io-error detected\n");
	}
	lengthRef=atoi(buf8);
	printf (" Original LengthRef %u , Dec %u length \n", lenRef, lengthRef);

	//3.-Length MapOP8 --> 8 bytes.. up to 99.999.999
/*	uint32_t lengthMap8O=0; //8 bytes
	buf8[8]={0};//=itoa(lengthReads,10);

	BZ2_bzRead(&bzerror, bfp, (char *)buf8, leng);

	if (bzerror == BZ_IO_ERROR){
	    	printf("bz-io-error detected\n");
		}
	lengthMap8O=atoi(buf8);
	printf (" Original LMap8 %u , Dec %u length \n", LenOpMap8, lengthMap8O);
*/
	//****4.-LEngth Bin Inst (posBInst) --> 8 bytes.. up to 99.999.999
	uint32_t lengthBI_O=0; //8 bytes
	buf8[8]={0};//=itoa(lengthReads,10);
	BZ2_bzRead(&bzerror, bfp, (char *)buf8, leng);

	if (bzerror == BZ_IO_ERROR){
	   	printf("bz-io-error detected\n");
	}
	lengthBI_O=atoi(buf8);
	printf (" Original LeBIN %u , Decomp %u length \n", posBInst, lengthBI_O);



	//********************************************READ READ LENGTH.
	char bufInt[5];
	uint32_t lengthReaD; //32 bits
	int len_I=5;

	BZ2_bzRead(&bzerror, bfp, (char *) bufInt, len_I);  //Variable de error, en el archivo que va, la data del buffer, longitud.
   	if (bzerror == BZ_IO_ERROR){
       	printf("bz-io-error detected\n");
   	}

	len_I=strlen(bufInt);
	lengthReaD=atoi(bufInt);
	printf (" Original %u , string %s length reaD %u\n", lengthReads, bufInt, lengthReaD);

	//*******************************************************Read MAP
	long Count=0;
	uint8_t buf2[BufCapacity+1];
	char *auxMap; //Reserve Memory
	uint8_t *MapDec=(uint8_t *)malloc((LenOpMap8+2)*sizeof(uint8_t));
	//redefine BufCapacity
	//auxMap= (char *)malloc((LenOpMap8+1)*sizeof(char)); //Pilas con este +1
	//auxMap[0]='\0';
	while ((bzerror == BZ_OK)&&(Count<= ((int)(LenOpMap8+1)-BufCapacity))){
	    BZ2_bzRead(&bzerror, bfp, (uint8_t *)buf2, BufCapacity);
	    							//BZ2_bzRead(&bzerror, bfp, (char*)buf1, BufCapacity);
	  	if (bzerror == BZ_IO_ERROR){        		printf("bz-io-error detected\n");    break;		}
	  								//  strcat(auxMap, buf1);//concatenate
	    memcpy(MapDec+Count, buf2, BufCapacity);//+1?

	  	if (DEBUG)for ( long g=Count;g<(Count+BufCapacity);g++){
	  			if (MAPAO_8[g]!=MapDec[g]) printf ("  %u %u  A_MISMATCH in Map at position %ld \n ",MAPAO_8[g],MapDec[g], g);
	  	}

	    Count=Count+BufCapacity;
	    //memset(buf1, 0, BufCapacity); //Re assign into final variable// Needed?
	}
	//memset(buf2, 0, BufCapacity); //Re assign into final variable// Needed?
	//******************************REST OF DATA
	diff = (LenOpMap8+1)-Count; 	    //procesar el resto de datos en el bbuffer
    //printf ("A Len %u Count %u diff %lu \n", LenOpMap8, Count, diff);
	if (diff>0) {
		BZ2_bzRead(&bzerror, bfp, (uint8_t *)buf2, diff);

		if (bzerror == BZ_IO_ERROR){	     			printf("bz-io-error detected\n");	  		 	}

		if (DEBUG) for ( long g=Count;g<(Count+diff);g++){//DUDA LenOpMap8 o (LenOpMap8+1)
				if (MAPAO_8[g]!=buf2[g-Count]) printf ("  %u %u  C_MISMATCH in Map at position %ld \n ",MAPAO_8[g],buf2[g-Count], g); //PILAS TODOS DAN MISMATCH
				//printf 4 posiciones de ambos vectores
		}

		memcpy(MapDec+Count, buf2, diff); //strcat(auxMap, buf1);//concatenate
	}
	//printf ("B Len %i Count %i diff %lu \n" ,LenOpMap8, Count, diff );
	//(compare to verify both maps)
	//Count, pilas, daña los primeros !
	if (DEBUG) for ( long g=Count;g<(LenOpMap8+1);g++){//DUDA LenOpMap8 o (LenOpMap8+1)
		if (MAPAO_8[g]!=MapDec[g]) printf ("  %u %u  BMISMATCH in Map at position %ld \n ",MAPAO_8[g],MapDec[g], g); //PILAS TODOS DAN MISMATCH
	} //O , LenOpMap8.

	//free (buf1);

	//ReadCantOfInstructions (Save before, it has not bee done)
	//************************************Read INSTRUCTIONS
	Count=0;
	uint8_t *DecBinInst=(uint8_t *)malloc((posBInst+1)*sizeof(uint8_t));
	//char *AuxDBI=(char *)malloc((posBInst+1)*sizeof(char)); //Reserve Memory
	//AuxDBI[0]='\0';;
	while ((bzerror == BZ_OK)&&(Count<= ((int)(posBInst+1)-BufCapacity))){
    	BZ2_bzRead(&bzerror, bfp, (uint8_t *)buf2, BufCapacity);
 		if (bzerror == BZ_IO_ERROR){        		printf("bz-io-error detected\n");    break;		}


 		memcpy(DecBinInst+Count, buf2, BufCapacity);

 		/*for ( long g=Count;g<(Count+BufCapacity);g++){
	  		if (BinInst[g]!=DecBinInst[g]) printf ("  %u %u  (BI)A_MISMATCH in Map at position %ld \n ",BinInst[g],DecBinInst[g], g);
	  	}*/

 		if ((DEBUG)&&(Count<6001)) for ( long g=Count;g<(Count+BufCapacity);g++){
 			if (BinInst[g]!=buf2[g-Count]) printf ("  %u %u  (BUF_BI)A_MISMATCH in Map at position %ld \n ",BinInst[g],buf2[g-Count], g);
 		}

		Count=Count+BufCapacity;

	}

	//memset(buf2, 0, BufCapacity);
	//leer remanente //******************************REST OF DATA
	diff = (posBInst+1)-Count; 	    //procesar el resto de datos en el bbuffer
	//printf ("psBInst %lu , Count %lu diff % lu", posBInst+1,Count, diff);
	if (diff>0){
		BZ2_bzRead(&bzerror, bfp, (uint8_t *)buf2, diff);
		if (bzerror == BZ_IO_ERROR){	     			printf("bz-io-error detected\n");	  		 	}

		memcpy(DecBinInst+Count, buf2, diff);

		if (DEBUG) for ( long g=Count;g<(Count+diff);g++){//DUDA LenOpMap8 o (LenOpMap8+1)
			if (BinInst[g]!=buf2[g-Count]) printf ("  %u %u  (Buf_BI)C_MISMATCH in Map at position %ld \n ",MAPAO_8[g],buf2[g-Count], g);
		}

	/*	for ( long g=Count;g<(Count+diff);g++){
			if (BinInst[g]!=DecBinInst[g]) printf ("  %u %u  (BI)A_MISMATCH in Map at position %ld \n ",BinInst[g],DecBinInst[g], g);
		}*/
	}


	//strcat(AuxDBI, buf1);
	//memset(buf1, 0, diff);
	//DecBinInst = (uint8_t *) AuxDBI; //Convert
	//(compare to verify both arrays)
	/*
	 * for (long g=0;g<(posBInst);g++){//DUDA posBInsto (posBInst+1)
		if (BinInst[g+1]!=DecBinInst[g]) printf ("MISMATCH in Map at position %ld", g); //PILAS TODOS DAN MISMATCH
	}
	 */


	BZ2_bzReadClose(&bzerror, bfp);
	if (fpin) fclose(fpin);
	//Free memory

	//ceil, tail
	printf ("LenRef %u, Lenref/8 %u , LenopMap8 %u \n", lenRef, lenRef/8,LenOpMap8 );



/************************************************END OF LOW LEVEL DECOMPRESSOR (pbzip2 HLI: High level interface (file output) )********************************************************************/


/************************************************ Start LOW LEVEL DECOMPRESSOR (pbzip2 LLI: LOW LEVEL INTERFACE)********************************************************************/


	/*   //allocate and initializa bz2_stream. bzallo, bzfree y opaque made null
	 *  //int BZ2_bzDecompressInit ( bz_stream *strm, int verbosity, int small );
	 *  //reserva de memoria para la data descomprimida ?
	//***len = total bytes to decompress? how will I know ?
	bytes_input=0; //ajgs
	while (bytes_input < len) {  //mientras haya data por descomprimir ?
	    isDone = false;

	    // Initialize the input buffer and its length
	    size_t in_buffer_size = len -bytes_input;
	    the_bz2_stream.avail_in = in_buffer_size;
	    the_bz2_stream.next_in = (char*)data +bytes_input; //***data is a pointer where the compressed data is read ?

	    size_t out_buffer_size =  output_size -bytes_uncompressed;  // size of output buffer //***Initialize intput bytes
	    if (out_buffer_size == 0) {  // out of space in the output buffer
	      break;
	    }

	    the_bz2_stream.avail_out = out_buffer_size;                    //
	    the_bz2_stream.next_out = (char*)output +bytes_uncompressed;  // ***output buffer, who is iteratively updated

	    ret = BZ2_bzDecompress(&the_bz2_stream); //acto de la descompression en si
	    if (ret == BZ_MEM_ERROR ) { 	      throw Bzip2Exception("Bzip2 Memory error ", ret); 	    }
	    if (ret != BZ_OK && ret != BZ_STREAM_END) { 	      throw Bzip2Exception("Bzip2 failed. ", ret); 	    }


        //debe cuidarse posibles errores de memoria al invocar  BZ2DEcompress. ver pag 16

	   bytes_input += in_buffer_size - the_bz2_stream.avail_in;
	   bytes_uncompressed += out_buffer_size - the_bz2_stream.avail_out;

	    *data_consumed =bytes_input;


	  }//end while bytes
*/


	//Bloque I (sequential)

	/*Hi level guide
	 *
	 *
	 */

	//1.- Read RefLocation= --> 32 bytes. (location). Dificultad: Manejar tamaño

	leng = 32;
	/*char*/ buf32[32]={0};

	//Que voy a comprimir, que tamaño tiene, donde lo descomprimo ?
	/*   //allocate and initialize bz2_stream. bzalloc, bzfree y opaque made null
		 *  //int BZ2_bzDecompressInit ( bz_stream *strm, int verbosity, int small );
		 *  //reserva de memoria para la data descomprimida ?
		//***len = total bytes to decompress? how will I know ?
		bytes_input=0; //ajgs
		    isDone = false;

		    // Initialize the input buffer and its length
		    size_t in_buffer_size = len -bytes_input;
		    the_bz2_stream.avail_in = in_buffer_size;
		    the_bz2_stream.next_in = (char*)data +bytes_input; //***data is a pointer where the compressed data is read ?

		    size_t out_buffer_size =  output_size -bytes_uncompressed;  // size of output buffer //***Initialize intput bytes
		    if (out_buffer_size == 0) {  // out of space in the output buffer
		      break;
		    }

		    the_bz2_stream.avail_out = out_buffer_size;                    //
		    the_bz2_stream.next_out = (char*)output +bytes_uncompressed;  // ***output buffer, who is iteratively updated

		    ret = BZ2_bzDecompress(&the_bz2_stream); //acto de la descompression en si
		    if (ret == BZ_MEM_ERROR ) { 	      throw Bzip2Exception("Bzip2 Memory error ", ret); 	    }
		    if (ret != BZ_OK && ret != BZ_STREAM_END) { 	      throw Bzip2Exception("Bzip2 failed. ", ret); 	    }


		   bytes_input += (in_buffer_size - the_bz2_stream.avail_in);
		   bytes_uncompressed += (out_buffer_size - the_bz2_stream.avail_out);

		    *data_consumed =bytes_input;

	*/


	//Bloque II (sequential)
		//Bloque III (sequential)
		//Bloque IV (sequential)
		//Bloque V (Iterative), por ajustar.

	/* **len = total bytes to decompress? how will I know ?
	bytes_input=0; //ajgs
	while (bytes_input < len) {  //mientras haya data por descomprimir ?
	    isDone = false;

	    // Initialize the input buffer and its length
	    size_t in_buffer_size = len -bytes_input;
	    the_bz2_stream.avail_in = in_buffer_size;
	    the_bz2_stream.next_in = (char*)data +bytes_input; //***data is a pointer where the compressed data is read ?

	    size_t out_buffer_size =  output_size -bytes_uncompressed;  // size of output buffer //***Initialize intput bytes
	    if (out_buffer_size == 0) {  // out of space in the output buffer
	      break;
	    }

	    the_bz2_stream.avail_out = out_buffer_size;                    //
	    the_bz2_stream.next_out = (char*)output +bytes_uncompressed;  // ***output buffer, who is iteratively updated

	    ret = BZ2_bzDecompress(&the_bz2_stream); //acto de la descompression en si
	    if (ret == BZ_MEM_ERROR ) { 	      throw Bzip2Exception("Bzip2 Memory error ", ret); 	    }
	    if (ret != BZ_OK && ret != BZ_STREAM_END) { 	      throw Bzip2Exception("Bzip2 failed. ", ret); 	    }


        //debe cuidarse posibles errores de memoria al invocar  BZ2DEcompress. ver pag 16

	   bytes_input += in_buffer_size - the_bz2_stream.avail_in;
	   bytes_uncompressed += out_buffer_size - the_bz2_stream.avail_out;

	    *data_consumed =bytes_input;


	  }//end while bytes
*/




		//Bloque VI (Iterative)
	     /*if (ret == BZ_STREAM_END) {                   //No tocar dejar como está
		      ret = BZ2_bzDecompressEnd(&the_bz2_stream);
		      if (ret != BZ_OK) { 	        throw Bzip2Exception("Bzip2 fail. ", ret); 	      }
		      isDone = true;
		    }*/


	/************************************************END OF LOW LEVEL DECOMPRESSOR (pbzip2 LLI)********************************************************************/




/***************************************************DECOMPRESSOR****************************************************************/

	//Reconstruct MAP
	//Ask for Memory
		uint8_t *MapRec;
		if (MapRec=(uint8_t *)calloc(lenRef+64,sizeof(uint8_t))) {printf (" Creado mapa ");}else printf (" No se puede crear mapa ");
		DeConstructOptimizedMap(MapRec , MAPAO, lenRef, LenOpMap);
		DeConstructOptimizedMap8(MapRec , MAPAO_8, lenRef, LenOpMap8);
		//free

		if (DEBUG) for (uint64_t g=0;g<=lenRef; g++){ //Just to verify
			if (MAPA[g]!=MapRec[g]) printf ("Mapa %i MAPA_RE %i g %ld Grave que mapas son diferentes \n",MAPA[g],MapRec[g],g);
		}

		if (DEBUG) if (map=fopen("MAPA_Reconst.txt","w")){
			for  (long y=1;y<=lenRef;y++) fprintf (map,"%i", MapRec[y]);
			fclose(map);
		}
		//if free

	//**************************************REFERENCE READER***********************************************/
			//char *ref, *PrevRead , *RefLocation; //UNCOMMENT ALL THIS
			//int lenRef; 	//ReferenceLength
			RefLocation="./SamData/lambda_virus.fa";
			ref= load_reference(RefLocation, &lenRef);
			//printf("longitud de la ref %i deberia ser 48502 y es %i \n", strlen(ref), lenRef-1);//En SAM la ref empieza en 1


	FILE *FInst=fopen("DecompressedInst.txt", "w"); //File for the decompressed instructions alone (not the final reads)

	uint32_t auxPosBInst=1;
	uint16_t *_Offsets, CountDels=0, CountIns=0,_CaEr;
	uint8_t  _MoreFrags,*_MoreErr,*_Oper,*_BaseRead,*_BaseRef, *DistReads,*_OPsCode, *_OPsBase, MapLocPrev;
	char _Match;				//Match==STRAND

	//Decompressor Output
	_Offsets=(uint16_t *)malloc(MAX_ERROR_READ*sizeof(uint16_t));  //Decompressed Offsets
	if (_Offsets==NULL) printf ("Not enough memory mfor Offsets %ld \n",i);

	_MoreErr=(uint8_t *)malloc(MAX_ERROR_READ*sizeof(uint8_t));	//Decompressed MoreErrors
	if (_MoreErr==NULL) printf ("Not enough memory for MoreError %ld \n",i);

	_BaseRead=(uint8_t *)malloc(MAX_ERROR_READ*sizeof(uint8_t));	//Decompressed Base Reads
	if (_BaseRead==NULL) printf ("Not enough memory for Base %ld \n",i);

	DistReads=(uint8_t *)malloc(MAX_ERROR_READ*sizeof(uint8_t));	//Auxiliar value to help in the developmenta stage to calculate distance
	if (DistReads==NULL) printf ("Not enough memory for DistReads %ld \n",i);

	_OPsCode=(uint8_t *)malloc(MAX_ERROR_READ*sizeof(uint8_t)); 	//Decompressed Operations
	if (_OPsCode==NULL) printf ("Not enough memory for Opcode %ld \n",i);

	uint64_t PrevPos, CurrentPos=0; //long AuxInd, AuxIndPrev; //uncomment when separate modules
	_MoreFrags=0;

	map=fopen("MAPA_GeneratedPositions.txt","w");

	for (long i=startMapped ; i<=tope; i++ ){//Only for mapped reads//Iterar desde StartMap hasta tope
			//En condiciones finales no se usa Indexpos nada sino directo la posicion origninal en q son despcomprimidos
			AuxInd=IndexPos[i];

		//*******************Decoding the binary instruction
		_BaseRef=BaseRef[AuxInd]; // Base References. During the development is used for verification purposes. //Decompressor Input only

		//MORALEJA TODO SE ESTA CALCULANDO BIEN; NO HAY ERROR EN CURRENTPOS NI EN POSITION NI EN AUXIND, QUIZAS EN MOREFRAGS
								//NEXT STEP ; COLOCAR TODO COMO DEBERIA SER EL SISTEMA DE INVOACIONES DEFINITIVO CON STARTMAP Y AUXIND y ENTONCES AJUSTAR
		PrevPos= CurrentPos;
		if ((_MoreFrags==0)) {
			CurrentPos=GetNextPos(PrevPos, MapRec, lenRef);
			if (DEBUG) fprintf (map,"%lu\n", CurrentPos); //This must be equal to the sorted original non repeated positions file
		}

/*Quit*/if (( Position[AuxInd]!=CurrentPos)) printf ("Original %i and decoded %lu positions are different in REad %lu \n",Position[AuxInd],CurrentPos,AuxInd);

		//DECODER
		if ((Match[AuxInd]!='U')&&(Match[AuxInd]!='u')&&(Position[AuxInd]!='0')){
			Bin2Inst(BinInst,  &auxPosBInst,  &_Match, &_MoreFrags, &_CaEr, _Offsets, _MoreErr,
		          _OPsCode, _BaseRead, _BaseRef, AuxInd, FInst, DistReads, &CountDels, &CountIns);

		//BUILDER
		int aux_;
		if (_Match =='f') _CaEr--; //PARA EL CASO FORWARD , dado que por la estructura del descriptor siempre arroja un error de mas /*cuidar si hay excepciones*/
			 ReadLength= longRead[AuxInd]; //No necesario en la version final.

			  aux_=BuildRead(_Match, CountDels, CountIns, ReadLength, CurrentPos, ref,
									   AuxInd, _Offsets, _CaEr,  _OPsCode, DistReads, Reads[AuxInd]);
		}

    }//// for

	fclose(map);

	if (_Offsets)    free (_Offsets);
	if (_MoreErr)     free (_MoreErr);
	if (_BaseRead)    free (_BaseRead);
	if (DistReads)	  free (DistReads);
	if (_OPsCode)	  free (_OPsCode);

	if (FInst) fclose (FInst);

	if (ref!=NULL) free(ref);
	if (MAPA) free (MAPA);
	if (MAPAO) free (MAPAO);
	//**Liberando la memoria usada para estructurar las instrucciones
	if (Position) free(Position);


	//free (Moerr, mediante un ciclo);
	/*for(int z=0; z<NUM_READSM;z++){ //YA liberado
		if (Offsets[z]) free(Offsets[z]);
		if (OPsCode[z]) free(OPsCode[z]);
		if (OPsBase[z]) free(OPsBase[z]);
	};*/


	/*	if ((Match[i-1]=='U')||(Match[i-1]=='u')){
					//strcpy(UnMapped[UC],Reads[i-1]);
					UnMapped[UC]=Reads[i-1];
					UC++;
				}*/


	end = clock();
	double time_spent1,time_spent2;

	time_spent1= (double)(end - begin) / CLOCKS_PER_SEC;
	printf("\n Tiempo TOTAL %f \n", time_spent1);

/***************************************************END OF DECOMPRESSOR****************************************************************/


	return(0);
}
/**************************************************MAIN END****************************************************************************/

/*******************************************************************************************************************************/
/******************************************************Reconstruction ( Instr to Read  )*****************************************/
/*******************************************************************************************************************************/


//Implementar el complementador de reads (T,A)(C,G)
void ComplementRead(char *Read, int length){
	char Compl;
	int i;
	for (i=0; i<length;i++){
		switch(Read[i]){
		    case 'A': Compl='T'; break;
		    case 'a': Compl='T'; break;
		    case 'C': Compl='G'; break;
		    case 'c': Compl='G'; break;
		    case 'G': Compl='C'; break;
		    case 'g': Compl='C'; break;
		    case 'T': Compl='A'; break;
		    case 't': Compl='A'; break;
		    case 'N': Compl='N'; break; //It is supposed not to find N's in the reference
		    case 'n': Compl='N' ;break;
		    default: printf ("**Error building the complement, wrong input base**");
		}

		Read[i]=Compl;
	}

};

//Implementar el Inversor de reads (T,A)(C,G)
void ReverseRead(char *Read, int length){
	char aux;
	for (int i=0; i<length;i++){
		aux=Read[length-i-1];
		Read[length-i-1]=Read[i];
		Read[i]=aux;
	}

};





/********************Reconstruccion de la Instruccion (Previamente probadas) **************/
/****This function builds and verify y verificar el read en funcion de las instrucciones generadas desde el decompresor*/
//En la version final debe quitarse y dejarse aparte la parte de verificacion.
//Depurar variables locales
//Obtener la long del read desde AFUERA (byte 0, 2 a la (byte+5)) //CUIDADO PORQUE EN ESTE MOMENTO NO ES ESTATICO y debe tomarse uno a uno!!!

//bloquear invocacion con (Distancias) en l read mayores a 4???  O al menos verificar?
//NOTA : SI EL CONSUMO DE TIEMPO ES EXAGERADO ESE MODULO SE PUEDE FUSIONAR CON EL QUE DESCOMPRIME
int BuildRead(char strand, uint16_t NDel, uint16_t NIns,  uint16_t LenRead, uint32_t MapPos, char *ref1, long NRead, uint16_t *Of,
		uint16_t lendesc,  uint8_t *Op, uint8_t *Distance,	char *Read /*Olg olg,	uint8_t NInst Instr *MInst, , int lenR, uint8_t Match,  uint8_t *Ba*/){
	//Nread numero identiiador del read
	uint16_t tam;
	char *auxR, *MiRead;

  	int w,iguales=0, posRead, posRef, ofAct,band=1, borradasF=0, borradasR=0, InsR=0, InsF=0; //Depurar las necesarias
  	static int raros=0;

  	FILE *aF, *aU, *DecB;

  	//1.Calcular la longitud del read a generar y guardarla en tam:

	tam=LenRead+NDel;

  	//CAlculo mi Read Inicial a partir de la referencia
  	auxR=ref1+MapPos;
  	if (MiRead=(char*)malloc(sizeof(char)*((tam)+1))){
  	 	MiRead[0]='\0';
  		strncat(MiRead, auxR, tam);
  	}else printf ("No hay memoria suficiente para crear MiRead");


  	////Aplicar transformaciones del caso REVERSO PERFECTO
  	////APlicar Transformaciones, Complemento  (sin reverso pues ya se incluye en el analisis del caso reverso.)
	//*EYE ON IT**if ((strand=='R')||(strand=='E')||(strand=='r')||(strand=='e'))  	     ReverseRead(MiRead, (int) tam);
  	if ((strand=='C')||(strand=='E')||(strand=='c')||(strand=='e'))  	     ComplementRead(MiRead, (int) tam);


  	//Solo usado para verificacion de ReadsRaros durante el desasrollo
  	if ((DEBUG)&&(DecB=fopen("DecompressedBuilder.txt","a"))){ //Weird Reads File
  		fprintf (DecB,"\n /******************INITIAL REF  READ %ld ****strand %c ******NInS %i *** NDELS %i ******************/", NRead, strand, NIns, NDel);
  		fprintf (DecB,"\n %s \n",MiRead);fflush(stdout);
  	}


  	//FACTORIZAR calculo de posicion y del AUXILIAR
  	if ((strand=='f')||(strand=='c')) { //forwardError, complementError(Always is Forward)
  				//Si es FORWARD:  evalúa el descriptor de izquierda a derecha, se ignora el ultimo entero del descriptor (el de mas a la izq) .  Se parean los números con el carácter de la izquierda. La primera instrucción en la ED se ejecuta de primera de derecha a izq sobre la REF.
  		char aux;
  		posRead=0;

  		char *RefB;
  		RefB=(char *)malloc(lendesc+1*sizeof(char));
  		//CalculateRefsF(RefB, MapPos, ref1, Of, lendesc,  Op); //**QUITAR***/

  		for (int y=0;y<=lendesc ;y++){//280317

 			switch(Op[y]) {
  			   case 's': // De acuerdo al MD sustitucion en en read, -->sustitucion en la referencia

  				   posRead=posRead+Of[y]; //Calculo la posicion del Offset en el read

  				   if (strand=='f') aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead+borradasF-InsF])+(int)Distance[y])%5);
  				   if (strand=='c') aux=DNAComplement(Index2Base(((int)BaseIndex(ref1[MapPos+posRead+borradasF-InsF])+(int)Distance[y])%5));

  				   if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB_Old %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);

				   MiRead[posRead]=aux;
				   posRead++;			  //Descuento la base actualizada

  			   break;
  			   case 'S': // De acuerdo al MD sustitucion en en read, -->Double contiguous substitution in the reference

  			   		posRead=posRead+Of[y]; //Calculo la posicion del Offset en el read

  			   		if (strand=='f') aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead+borradasF-InsF])+(int)Distance[y])%5);
  			   	  	if (strand=='c') aux=DNAComplement(Index2Base(((int)BaseIndex(ref1[MapPos+posRead+borradasF-InsF])+(int)Distance[y])%5));

  			   	    if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB_Old %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);
  			   		MiRead[posRead]=aux;
  			   		posRead++;			  //Descuento la base actualizada

  			 		MiRead[posRead]=aux;   //Contiguous base, Does not need to recalculate the aux
  			   		posRead++;			  //Descuento la base actualizada

  			   	break;

  			   case 'd': // De acuerdo al MD delecion en la referencia -->insercion en el read al final
  			  		posRead=posRead+Of[y]; //Calculo la posicion del Offset en el read
  			  		aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead+borradasF-InsF])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia . (REVISAR)
  			  	    if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y], '-', aux, Distance[y]);

  					for(int q=posRead;q<tam-1;q++){
  						  MiRead[q]=MiRead[q+1];
  					}
  					borradasF++;

  			   break;
  			   case 'D': // De acuerdo al MD delecion en la referencia -->insercion en el read al final
  			   		posRead=posRead+Of[y]; //Calculo la posicion del Offset en el read
  			   		aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead+borradasF-InsF])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia . (REVISAR)
  			   	    if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y], '-', aux, Distance[y]);
  			   		for(int q=posRead;q<tam-2;q++){
  			   			  MiRead[q]=MiRead[q+2];
  			   		}
  			   		borradasF=borradasF+2;

  			   break;

  			 case 'T': // De acuerdo al MD delecion en la referencia -->insercion en el read al final
  			   		posRead=posRead+Of[y]; //Calculo la posicion del Offset en el read
  			   		aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead+borradasF-InsF])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia . (REVISAR)
  			   	    if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y], '-', aux, Distance[y]);
  			   		for(int q=posRead;q<tam-3;q++){
  			   			  MiRead[q]=MiRead[q+3];
  			   		}
  			   		borradasF=borradasF+3;

  			  break;

  			 case 'C': // De acuerdo al MD delecion en la referencia -->insercion en el read al final
  			  		posRead=posRead+Of[y]; //Calculo la posicion del Offset en el read
  			  		aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead+borradasF-InsF])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia . (REVISAR)
  			  	    if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y], '-', aux, Distance[y]);
  			   		for(int q=posRead;q<tam-4;q++){
  			   			  MiRead[q]=MiRead[q+4];
  			   		}
  			   		borradasF=borradasF+4;

  			  break;

  			   case 'i': // De acuerdo al MD una insercion -->delecion
  				   //Inserte en el read para que se asemeje a la referncia, corresponde a eliminar en la referencia (MiRead)
  				   //la insercion ocurre en l a POS SIGUIENTE, no en la exacta

  				   posRead=posRead+Of[y]; //Calculo la posicion del Offset en el read


  				   if (strand=='f') aux=Index2Base((int)Distance[y]);
  				   if (strand=='c') aux=DNAComplement(Index2Base((int)Distance[y]));


  				   if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y], '-', aux, Distance[y]);

  				   for(int q=(tam-1);q>posRead;q--){
  					  MiRead[q]=MiRead[q-1];
  				   }

  				   MiRead[posRead]=aux;
  				   posRead++;
  				 InsF++;
  	 			break;
  			    case 'I': // De acuerdo al MD una double insertion -->Double deletion in the REF
  			   				   //Inserte en el read para que se asemeje a la referncia, corresponde a eliminar en la referencia (MiRead)
  			   				   //la insercion ocurre en l a POS SIGUIENTE, no en la exacta

  			   		posRead=posRead+Of[y]; //Calculo la posicion del Offset en el read
  			   		for(int q=(tam-1);q>posRead;q--){
  			   			  MiRead[q]=MiRead[q-1];
  			   		}
  			   	    if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y], '-', 'N', Distance[y]);
  			   		MiRead[posRead]='N';
  			   		posRead++;

  			   	InsF++;
  			   	break;
  				case '_': case '?': case 'X':
  				break;
  			}//end switchF

  		}//endForyF
  		free(RefB);
  	}else{  	//REVERSE and Reverse Complement,
  		if ((strand=='r')||(strand=='e')) {
  		//Si es Reverse: Se evalúa el descriptor de derecha a izquierda, se ignora el primer entero del descriptor (el de mas a la izq) .  Se parean los números con el carácter de la izquierda. La primera instrucción en la ED se ejecuta de primera de derecha a izq sobre la REF.
          posRead=(tam-1)-NIns;//+Ndel;

        char aux;//the Nucleotide Base
  		for (int y=0;y<lendesc;y++){

  			switch(Op[y]) {
  			   case 's': // De acuerdo al MD sustitucion en read, -->sustitucion en la referencia

  				   posRead=posRead-Of[y];//+borradasR //Calculo la posicion del Offset en el read

  				   if (strand=='r') aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia .
  				   if (strand=='e') aux=DNAComplement(Index2Base(((int)BaseIndex(ref1[MapPos+posRead])+(int)Distance[y])%5));

  				   if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);
  				   MiRead[posRead]=aux;
  				   posRead--;			  //Descuento la base actualizada

  			   break;
  			   case 'S': // De acuerdo al MD Double sustitucion en read, -->sustitucion en la referencia

  				   posRead=posRead-Of[y]; //+borradasR//Calculo la posicion del Offset en el read

  				   if (strand=='r') aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia .
  				   if (strand=='e') aux=DNAComplement(Index2Base(((int)BaseIndex(ref1[MapPos+posRead])+(int)Distance[y])%5));


  				   if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);
 				   MiRead[posRead]=aux;
  				   posRead--;			  //Descuento la base actualizada

  				   //Contiguous BASE
  				   MiRead[posRead]=aux;
  				   posRead--;			  //Descuento la base actualizada

  			   break;
  			   case 'd': // De acuerdo al MD delecion en la referencia -->insercion en el read
  				   posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read
  				   aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia . (REVISAR)
  				   if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);
    			   for(int q=posRead;q<tam-1;q++){ /***********ver si este tope con strlen o incluir NINS************/
  					    MiRead[q]=MiRead[q+1];
  				   }

  				  posRead--;
  				  borradasR++;
  			   break;
  			   case 'D': // De acuerdo al MD Double delecion en la referencia -->insercion en el read
  				   posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read
  				   aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia . (REVISAR)

  				   if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);
    			   for(int q=posRead-1;q<tam-2;q++){ /***********ver si este tope con strlen o incluir NINS************/
  					    MiRead[q]=MiRead[q+2];
  				   }

  				  posRead=posRead-2;
  				  borradasR=borradasR+2;
  			   break;
  			   case 'T': // De acuerdo al MD Triple delecion en la referencia -->insercion en el read
  			  		posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read
  			  	    aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia . (REVISAR)

  			  	    if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);
  			    	for(int q=posRead-2;q<tam-3;q++){ /***********ver si este tope con strlen o incluir NINS************/
  			  		    MiRead[q]=MiRead[q+3];
  			  		}

  			  		 posRead=posRead-3;
  			  		 borradasR=borradasR+3;
  			   break;

  			   case 'C': // De acuerdo al MD Quadruple delecion en la referencia -->insercion en el read
  			  		posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read
  			  		aux=Index2Base(((int)BaseIndex(ref1[MapPos+posRead])+(int)Distance[y])%5);//SUma pos de la Ref+Distancia . (REVISAR)

  			  	    if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);
  			    	 for(int q=posRead-3;q<tam-4;q++){ /***********ver si este tope con strlen o incluir NINS************/
  			  		    MiRead[q]=MiRead[q+4];
  			  		 }

  			  		 posRead=posRead-4;
  			  		 borradasR=borradasR+4;
  			   break;

  			   case 'i': // De acuerdo al MD una insercion -->delecion
  				   //Inserte en el read para que se asemeje a la referncia, corresponde a eliminar en la referencia (MiRead)
  				   //la insercion ocurre en l a POS SIGUIENTE, no en la exacta

  				   posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read

  				   if (strand=='r') aux=Index2Base((int)Distance[y]);
  				   if (strand=='e') aux=DNAComplement(Index2Base((int)Distance[y]));

  				   if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, Distance[y]);

  				   for(int q=(tam-1);q>posRead+1;q--){   //utilizar la longitud de la variable mas quela del read
  					  MiRead[q]=MiRead[q-1];
  				   }
  				   MiRead[posRead+1]=aux;  //.Puede que toque sumar 1
  				   InsR++;

  			   break;

  			   case 'I': // De acuerdo al MD una Double insercion -->delecion
  			   				   //Inserte en el read para que se asemeje a la referncia, corresponde a eliminar en la referencia (MiRead)
  			   				   //la insercion ocurre en l a POS SIGUIENTE, no en la exacta

  			   	   posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read
  			   	   if (DEBUG) fprintf(DecB, "#Error %i Offset %u  Oper %c  RefB %c ReadBase %c Dist %i  \n", y+1,Of[y], Op[y],   ref1[MapPos+posRead], aux, 'N');
  			   	   for(int q=(tam-1);q>posRead+1;q--){   //utilizar la longitud de la variable mas quela del read
  			   		   MiRead[q]=MiRead[q-1];
  			   	   }
  			   	   MiRead[posRead+1]='N';

   				   InsR++;

  			   break;

  			   case '_': case '?': case 'X': case 'x'://posRead=posRead-Of[y];//Solo informacion de offset
  			   break;

  			}//end switch

  		}//endFory
  		}//end if (Match== r )||(Match==e)
  	}//Endelse Reverse
  	MiRead[LenRead]='\0';
  	fprintf(DecB, " \n");
  	//Verificacion de la reconstruccion
  	if (((DEBUG))&&(strcmp(Read,MiRead)!=0)){


  			fprintf (DecB,"/Original**********************WEIRD PAR FINAL ********lendesc %i *********************OlG, MiRead/ \n", lendesc);
  			fprintf (DecB,"\n Caso raro de error en sentido %c Builder read numero : %ld \n",strand, NRead);

  			fprintf (DecB,"\n %s\n",Read);
  			fprintf (DecB,"\n %s\n",MiRead);
  			raros++;
  			printf("\n BUILDER RAROS %i, read %ld, NiNs %i, Ndels %i y posRead %i editD %i ",  raros, NRead,   NIns, NDel, posRead, lendesc);
  	}
  	fclose(DecB);
  	return (strcmp(Read,MiRead));//Compara el READ generado con el del olg
};

//DNA COMPLEMENTS
char DNAComplement(char Base){
	switch (Base){
		case 'A': case 'a': return ('T'); break;
		case 'C': case 'c': return ('G'); break;
		case 'G': case 'g': return ('C'); break;
		case 'T': case 't': return ('A'); break;
		case 'N': case 'n': return ('N'); break;
	}

}


//Get Mapping position
uint64_t GetNextPos(uint64_t PrevPos, uint8_t *Map, long lenRef){
   PrevPos++;
   while ((PrevPos<lenRef)&& (Map[PrevPos]==0)){
	   PrevPos++;
    }
	return(PrevPos);
 };

/*******************************************************************************************************************************/
/******************************************************End of Reconstruction ( Instr to Read  )************************************/
/*******************************************************************************************************************************/

/*******************************************************************************************************************************/
/********************************************** Decompression ( Binary to Instr )***********************************************/
/********************************************(Convert Binary into Instructions) ************************************************/
/*******************************************************************************************************************************/

//Counts the number of bases corresponding to a Insertion oper (single/double)
int CountInsBases(uint8_t Oper){
	if (Oper=='i') return (1); else return (2);

}
//Counts the number of bases corresponding to a Deletion oper (single/double/triple/quad)
int CountDelBases(uint8_t Oper){
	switch (Oper){
		case 'd': return (1); break;
		case 'D': return (2); break;
		case 'T': return (3); break;
		case 'C': return (4); break;
	}

}

//Pendiente: Inicializar adecuadamente posBInst++ cuando este modulo se desincronize de la compresion.
	//cuidar y verificar sentido de la reconstruccion de la instruccion
//Input: An array (BinInst) where all the bytes were previously stored, and a position to know where this read begins
//This module is invoked for each read (each instruction), If the read is not unmapped.
//Output: Convertir Binario a  instruccion, recibe toda la informacion de un mismo read y genera la correspondiente codificacion binaria
//por la forma en que esta programado el compresor , el descompresor la retornara en el orden correcto de lectura, esto es
//los de sentido R en su orden natural, y los de F en el inverso, con los Offsets correctamente colocados en ambos casos

void Bin2Inst(uint8_t *BinInst,  uint32_t *posBInst,  char *strand, uint8_t *MoreFrags, uint16_t *lendesc, uint16_t *_Offsets,
		uint8_t *MoreErr, uint8_t *_Oper, uint8_t *BaseRead, uint8_t *BaseRef, uint32_t Nread, FILE *FInst,
		uint8_t *DistReads, uint16_t *CountDels, uint16_t *CountIns){ 		//Nread , DistReads and BaseRef are here only for the test phase
	char auxStrand;
	uint8_t rest=0x0, aux=0,MoreEr, Operation, ReadBase, auxMoreFrags;
	uint32_t auxPosInst=*posBInst;//+1; //Este mas 1 es necesario durante el desarroollo, en la version real de la descompresion debe autoregularse el punto exacto
	*CountDels=0; *CountIns=0;

	Bin2Preambulo(BinInst[auxPosInst],&auxMoreFrags, &auxStrand);//pilas paramtros //Meter en los respectivos arreglos que corresponden.
	*strand=auxStrand; *MoreFrags=auxMoreFrags;
		//printf(" read %i STRND %c MoreFrags %i POSICION %ld valor %i\n", Nread, *strand, *MoreFrags, auxPosInst, BinInst[auxPosInst]);
	if (DEBUG) fprintf(FInst, "NumRead %i MoreFrags %i STRAND %c \n", Nread, *MoreFrags,*strand );//PRINT2 FILE
	auxPosInst++;


		if (((*strand=='f')||(*strand=='r')||(*strand=='c')||(*strand=='e'))){ //If there are errors
		uint8_t  x=0; //To Count the errors
		MoreEr=1;  //Since all the letters are in lowercase, we know that there is at least one more error

		if ((*strand=='e')||(*strand=='r')){ 		//All reverse cases
			while (MoreEr==1){ //Each contiguous error is processed
					//printf("Before Offsets R read %i error  %i posicion %ld \n", Nread,x, auxPosInst); fflush (stdout);
				_Offsets[x]= Bin2Offset(BinInst,auxPosInst);
				if (DEBUG) fprintf(FInst, "#Error %i Offset %u     ", x+1,_Offsets[x] );
				auxPosInst++;

					//printf("Before 3rd R read %i error  %i posicion %ld valor %i\n", Nread,x, auxPosInst,BinInst[auxPosInst]); fflush (stdout);
				Byte3Bin2Inst(BinInst[auxPosInst],&MoreEr, &Operation, &ReadBase);//ReadBase is basically the distance
					//printf (" MoreEr %i Operation %i ReadBase %i \n",MoreEr, Operation, ReadBase);
				MoreErr[x]=MoreEr; _Oper[x]=Operation; DistReads[x]=ReadBase;

				if ((Operation=='i')||(Operation=='I')){
					*CountIns=*CountIns+1;//CountInsBases(Operation);
				}else{
					if ((Operation=='d')||(Operation=='D')||(Operation=='T')||(Operation=='C')){
						*CountDels=*CountDels+CountDelBases(Operation);
					}
				}; //printf ("_Oper[x] %i CountDels %i  CountIns %i ",_Oper[x], *CountDels, *CountIns);

				if ((Operation!='i')&&(Operation!='I')){
					//se puede eliminar el calculo siguiente para casos de delecion tambien
					BaseRead[x]=uint8_t (((int)BaseRef[x]+(int)DistReads[x])%5);      //Both bases must come as uint8
				}else{
					if (Operation=='i') BaseRead[x]=(int)DistReads[x];
						else BaseRead[x]='N';
				}

				if (DEBUG) fprintf(FInst, "MoreError %u Oper %c DistBas %u, BaseRef %c, BaseRead %c \n", MoreErr[x], Myitoa(_Oper[x]), DistReads[x],Index2Base(BaseRef[x]) ,Index2Base(BaseRead[x]) );
				auxPosInst++;

				x++;

			}//EndWh
		}else{  //All Forward cases
			while (MoreEr==1){ //Each contiguous error is processed
					//printf("Before Offsets F read %i error  %i posicion %ld \n", Nread,x, auxPosInst); fflush (stdout);
				_Offsets[x]= Bin2Offset(BinInst,auxPosInst);
				if (DEBUG) fprintf(FInst, "#Error %i Offset %u     ", x+1,_Offsets[x] );
				auxPosInst++;

					//printf("Before 3rd R read %i error  %i posicion %ld valor %i\n", Nread,x, auxPosInst,BinInst[auxPosInst]); fflush (stdout);
				Byte3Bin2Inst(BinInst[auxPosInst],&MoreEr, &Operation, &ReadBase);//ReadBase is basically the distance
					//printf (" MoreEr %i Operation %i ReadBase %i \n",MoreEr, Operation, ReadBase);
				MoreErr[x]=MoreEr; _Oper[x]=Operation; DistReads[x]=ReadBase;

				if ((Operation=='i')||(Operation=='I')){
									*CountIns=*CountIns+1;//CountInsBases(Operation);
				}else{
					if ((Operation=='d')||(Operation=='D')||(Operation=='T')||(Operation=='C')){
						*CountDels=*CountDels+CountDelBases(Operation);
					}
				}; //printf ("_Oper[x] %i CountDels %i  CountIns %i ",_Oper[x], *CountDels, *CountIns);

				if ((Operation!='i')&&(Operation!='I')){
				//se puede eliminar el calculo siguiente para casos de delecion tambien
					BaseRead[x]=uint8_t (((int)BaseRef[x]+(int)DistReads[x])%5);      //Both bases must come as uint8
				}else{
					if (Operation=='i') BaseRead[x]=(int)DistReads[x];
						else BaseRead[x]=0x4;//case for N
				}

				if (DEBUG) fprintf(FInst, "MoreError %u Oper %c DistBas %u, BaseRef %c, BaseRead %c \n", MoreErr[x], Myitoa(_Oper[x]), DistReads[x],Index2Base(BaseRef[x]) ,Index2Base(BaseRead[x]) );
				auxPosInst++;

				x++;

			}//endwhile
		}//Endelse
		*lendesc=x;  //PILAS ESTA ARROJANDO SEGMENTATION FAULT, EN CUAL CASO ???, read 32, F
	}//If Errors
	*posBInst=auxPosInst; //NOTA esta instruccion debe IR descomentada en la Version final del descompresor para que la posicion se autoregule
							//en esta version intermedia viene dada por la posicicion de compresion
	fprintf(FInst, "\n");
};//EndBin2Inst


//Convierte el primer byte de la cadena binaria en los Datos del preambulo.
void Bin2Preambulo(uint8_t FstByte, uint8_t *moreFrags, char *strand){
	uint8_t mask=0x8, aux=0x7; //1000, 0111 RESPECTIVELY

	//1. Considerar caso de dos preambulos perfectos contiguos. si esto se implementa?,

	if ((mask&FstByte)>0) { //MoreFrags *RESOLVER
		*moreFrags=1;
	} else *moreFrags=0;

	aux=aux&FstByte;

	switch(aux){   	 //Calculando Matches / STRAND, in the 3 LSB
		case 0x0: *strand='F';break;//Mayusculas los casos de matching PERFECTOS
		case 0x1: *strand='R';break;
		case 0x2: *strand='C';break; //010
		case 0x3: *strand='E';break; //011
		case 0x4: *strand='f';break;//Minusculas los casos de matching con Error//100
		case 0x5: *strand='r';break;//101
		case 0x6: *strand='c';break;//110
		case 0x7: *strand='e';break;//111
	}
	//printf ("Fbyte %i aux %i MF %i ST %c ",FstByte, aux, *moreFrags, *strand);
};

//Calculates the offset using the 2nd byte and the 2 MSB of the 3rd byte
uint16_t Bin2Offset(uint8_t *BinInst, uint32_t PosInst){
	uint8_t mask;
	uint16_t aux;
	mask=0xC0;// //0x11000000;
	aux=(mask&BinInst[PosInst+1]); //2MSB del 3er byte
	aux=aux<<2;
	return(aux|BinInst[PosInst]);
};

//Toma el tercer byte y lo ocnvierte en la respectiva data de instruccion
//2bits mas significativos: 2 bits mas significativos del offset cuando este es >=256
//3er bit mas significativo: Bit del MoreError
//4to, 5to y 6to bit mas significativo: OPER
//2 bits menos significativos: Base

void Byte3Bin2Inst(uint8_t Byte3,uint8_t *MoreErr, uint8_t *oper, uint8_t *DistReads){
	uint8_t mask,aux ;

	//Byte3=Byte3<<2; //DEscarto los dos bits del rest
	mask=0x20;//100000 //Creo la Mascara de More Error.     //CalcularMoreError

	if ((mask&Byte3)==mask) *MoreErr=1; else *MoreErr=0;
   // Byte3=Byte3<<1;

	mask=0x1c;//00011100;						//Calcular Oper
	aux=mask&Byte3;
	aux=aux>>2;

	switch(aux){
		case 0x0: *oper='s' ; break; // Sust Simple
		case 0x1: *oper='i' ; break; // InsSimple
 		case 0x2: *oper='d' ; break; // Delecion Simp  01000
		case 0x3: *oper='D' ; break; // Delecion Doble 01100
		case 0x4: *oper='T' ; break; // Delecion Triple 10000
		case 0x5: *oper='C' ; break; // Delecion cuadruple 10100
		case 0x6: *oper='S' ; break; // Sust Doble repetida contigua 11000
		case 0x7: *oper='I' ; break; // en N
	}

	mask=0x3;						//Calcular Base
	aux=mask&Byte3;

	if (*oper!='i') *DistReads=aux+1; //+1, +2, +3,+4 , para luego aplicar el mod 5, por el momento no se puede saber sin la ref
		else *DistReads=aux;
};
char Myitoa(uint8_t val){
char oper;
	switch(val){
			case 115: oper='s' ; break; // Sust Simple
			case 105: oper='i' ; break; // InsSimple
	 		case 100: oper='d' ; break; // Delecion Simp  01000
			case 68: oper='D' ; break; // Delecion Doble 01100
			case 84: oper='T' ; break; // Delecion Triple 10000
			case 67: oper='C' ; break; // Delecion cuadruple 10100
			case 83: oper='S' ; break; // Sust Doble repetida contigua 11000
			case 73: oper='I' ; break; // Ins en N
			default: oper='-' ;
		}
	return (oper);

}


/*******************************************************************************************************************************/
/**********************************************END of Decompression ( Binary to Instr  )*****************************************/
/*******************************************************************************************************************************/


//this prints a byte to the output
void PrintBinary2(unsigned int num, int nbits){
    for( int i = nbits-1 ; i>=0 ; --i )
        if( num & (1 << i) ) printf("1");
        else printf("0");
}
//this prints a byte to a file
void PrintByte2File(FILE *outB, uint8_t num, int nbits){
	for( int i = nbits-1 ; i>=0 ; --i )
	    if( num & (1 << i) ) fprintf(outB,"%c", '1');
	        else fprintf(outB,"%c", '0');
	 fprintf(outB,"\t");
};

/*******************************************************************************************************************************/
/******************************************************Compression ( Instr to Binary  )*****************************************/
/*******************************************************************************************************************************/

//Crea una mascara de Long bits de Der a Izq,  a partir de la posicion menos significativa , MSB[0][000000][1] LSB
uint8_t CreateMask8B(uint8_t pos_, uint8_t Long){
   uint8_t mask=0x01, aux=0x0, i;
   int pos=pos_;
   if (((pos<=8)&&(pos-(int)Long))>=0){
	   mask=mask<<(pos-1); aux=aux|mask;
		for(i=1; i<Long; i++){
			mask=mask>>1;  //aux=aux<<1;//<<
			aux=aux|mask;
		}
		if (pos-(int)Long>0) aux=aux<<(pos-(int)Long); //Introducir los bits restantes
   }
   return (aux);
}

//OJO Inicializar posBInst++;, se invoca este modulo para cada olg, es decir Instuccion
//if el read no es unmapped, y si tengo espacio suficiente een la estructura BinInst (Cuanto? espacio)

//This converts an instruction to binary code [Prelude] ([Offset]|Operation)
//Input: recibe toda la informacion de un mismo read y genera la correspondiente codificacion binaria
//Output: PLEASE BE CAREFUL ABOUT THE  DIRECTION used to return the Forward instructions, (RIght to left)or (LtoR (the opposite))
void Inst2Bin(uint8_t *BinInst,  uint32_t *posBInst, char strand, uint8_t MoreFrags, uint16_t lendesc, uint16_t *Offsets,
		uint8_t *Oper, uint8_t *BaseRead, uint8_t *BaseRef, FILE *outB, FILE *outb2c, FILE *bin, long i, uint16_t EdDis){
	
	uint32_t auxPosInst=*posBInst ;
	uint8_t rest=0x0, aux=0, MoreErr=1;
	int aux_i;

	struct rec my_record;//AUX for storing binary files

	auxPosInst++;
	BinInst[auxPosInst]=Preambulo(MoreFrags, strand , lendesc, EdDis);


	/*  //To verify correctness of the byte string for the preamble
	FILE *PAmbFile=fopen("PreambFile.txt", "a");
	fprintf(PAmbFile,"Readnumber %ld , posBInst %ld , byte %u ,  MORE FRAGS %i STRAND %c lendesc %i\n\n",i, auxPosInst,BinInst[auxPosInst],MoreFrags, strand , lendesc);
	if (PAmbFile) fclose(PAmbFile);
	 */

	PrintByte2File(outB, BinInst[auxPosInst], 8); //Print to a Bytes File (as a string of '0' and '1' )
	if (DEBUG) fprintf(outb2c,"%u ",BinInst[auxPosInst]);	  //Print to Unsigned Int File

	//char buffer;                                //Trying to print the binary value as a simple ASCII caracter, many ways tried, they did not work
	//itoa(BinInst[auxPosInst],&buffer,10); 	  //printf("%c", buffer);

	my_record.x=BinInst[auxPosInst];				//used to build the binary file
	if (DEBUG) fwrite(&my_record, sizeof(struct rec), 1, bin); //Print to bin FIle

	//If this is not a perfect matching
	//NOTA IMPORTANTE, EN REALIDAD AQUI SOLO VAN r o f en MINUSCULAS, pero por fase de PRUEBAS se colocara en mayusculas in order to satisfy the condition by the moment
    if ( ((lendesc>0)&&((strand=='r')||(strand=='e'))) || ((lendesc>1)&&((strand=='f')||(strand=='c'))) ){ //Si hay errores
		if ((strand=='r')||(strand=='e')){

			for (uint8_t  u=0;u<lendesc; u++){ 			//Converting each separated error of the read
				if ((Oper[u]!='_')&&(BaseRead[u]>=0)&&(BaseRead[u]<=4)){
					auxPosInst++;

					BinInst[auxPosInst]= Offset(Offsets[u], &rest);
					PrintByte2File(outB, BinInst[auxPosInst], 8);
					if (DEBUG) fprintf(outb2c,"%u ",BinInst[auxPosInst]);     //unsigned int file
					my_record.x=BinInst[auxPosInst];//for input and outputtests
					if (DEBUG) fwrite(&my_record, sizeof(struct rec), 1, bin);//binary file
out
					auxPosInst++;
	//if (((BaseRead[u]=='n')||(BaseRead[u]=='N'))&&(Oper[u]!='i')) printf ("****AAAAA     Ins en N ****");
					BinInst[auxPosInst]= TrdBitInst(u, rest, Oper, BaseRead, BaseRef[u], Offsets, lendesc, strand, &aux_i);
					PrintByte2File(outB, BinInst[auxPosInst], 8);
					if (DEBUG) fprintf(outb2c,"%u ",BinInst[auxPosInst]);
					my_record.x=BinInst[auxPosInst];
					if (DEBUG) fwrite(&my_record, sizeof(struct rec), 1, bin);//binary file

					u=aux_i;
				}//endif Oper
			}//endfo
		}else{  //Forward cases
			/*(in the Forward cases , the instructions must be processed from)
			 * Se evalúa el descriptor de izq a derecha se ignora el ultimo entero  entero (el de mas a la derecha).
			 *Se parean los números con el carácter de la derecha. La ultima instrucción de la ED se ejecuta de primera (de izq a der) sobre la ref.
			 * */

/***********///cuidar operaciones contiguas, salto de la referencia y de todo en operaciones de contiguidad, valor de mor errors por el sentido
				//Intentar guardarlas exactamente como se generan pero al contrario . podria ser con lendesc-1
			for (int  u=lendesc-1;u>=0; u--){ 		//Conversion de cada uno de los erroes (***lendesc o lendes-1)
				//SALTAR la base de la ref cuando aplica la contigua
//NO CHANGES YET
				if ((Oper[u]!='_')&&(BaseRead[u]>=0)&&(BaseRead[u]<=4)){
					auxPosInst++; //El i+1 solo se aplica en los casos de OFFSETS OK

					BinInst[auxPosInst]= Offset(Offsets[u+1], &rest); //[aqui se rueda el desplazamiento al elemento de la DERECHA por eso empiezo al reves]
					/*fprintf(outB,"B2f ");*/PrintByte2File(outB, BinInst[auxPosInst], 8);
					if (DEBUG) fprintf(outb2c,"%u ",BinInst[auxPosInst]);
					my_record.x=BinInst[auxPosInst];
					if (DEBUG) fwrite(&my_record, sizeof(struct rec), 1, bin);//bin only file

					auxPosInst++;
					//printf ("CALCULADO BYTE 2 F PosVec %i \n",auxPosInst);fflush(stdout);
					//if ((BaseRead[u]==BaseRef[u])&&(BaseRead[u]!='-')) printf ("***Read %ld F BASE ERROR %i ***\n", i, u+1);
					//fprintf(outB,"\n Oper %c BRead %u BRef %u lendesc %u VarI %i \n",Oper[u], BaseRead[u],BaseRef[u], lendesc,u);
/*ojo lendesc-1*/	BinInst[auxPosInst]= TrdBitInst(u, rest, Oper, BaseRead, BaseRef[u], Offsets, lendesc, strand,&aux_i);
					/*fprintf(outB,"B3f ");*/PrintByte2File(outB, BinInst[auxPosInst], 8);
					if (DEBUG) fprintf(outb2c,"%u ",BinInst[auxPosInst]);
					my_record.x=BinInst[auxPosInst];
					if (DEBUG) fwrite(&my_record, sizeof(struct rec), 1, bin);//bin only file

					u=aux_i;
					//printf ("CALCULADO BYTE 3 F Posvec %i  valor de u %i\n",auxPosInst, u);fflush(stdout);
				}//endif
			}//endfor


		}//ende_else
    }//endIF
	*posBInst=auxPosInst;
	fprintf(outB,"\n \n");
};//EndInst2Bin

//2 MSB (aka Most significant bit ): Offset's 2 MSB bits when this is >=256
//3er bit mas significativo: Bit del MoreError
//4to, 5to y 6to bit mas significativo: OPER
//2 bits menos significativos: Base
uint8_t TrdBitInst( int i, uint8_t  rest, uint8_t  *Oper, uint8_t  *BaseRead, uint8_t BaseRef,
		uint16_t *offset, uint16_t lendesc , char strand, int *aux_i){
	uint8_t mask=0x01, aux=0x0, mask2=0x0;
	int ii=i;
	aux=aux|rest;						  ////CalcularOffset,Adiciono los 2 ultimos bits del Resto anterior
    aux=aux<<1;

  //  printf("ii %i i %i auxi %i\n",ii, i, *aux_i);

    if ((strand=='R')||(strand=='e')) mask2=BitsOperR(Oper, BaseRead, offset, lendesc, &i); /*Returns the oper and the Redbase binary values*/
    if ((strand=='F')||(strand=='c')) mask2=BitsOperF(Oper, BaseRead, offset,  &i);

    //NOTA IMPORTANTE, EN REALIDAD AQUI SOLO VAN r o f en MINUSCULAS, pero por fase de PRUEBAS se colocaran en mayusculas
       	//para que se satisfaga la condicion
   	if (((strand=='R')||(strand=='e'))&&((i<lendesc-1)&&(Oper[i+1]!='_'))) aux=mask|aux;  ////Calculates MorErrors,
   		else if (((strand=='F')||(strand=='c'))&&(i>0)) aux=mask|aux;

   	//if ((strand=='R')&&(lendesc==5)) printf("AUX %u, oper %u \n", aux, mask2);

   	aux=aux<<3;
    aux=aux|mask2;
	aux=aux<<2;

	if((mask2==0x2)||(mask2==0x3)||(mask2==0x4)||(mask2==0x5)||(mask2==0x7)) {//esDelecion en la REF, o ins N
		mask=0x00;									//The base in the read does not really matters
	}else {
		if ((mask2!=0x1)){ //Is a single or double subst
			mask=BitsBase(BaseRead[ii], BaseRef);				//CAlculates the ReadBase
		}else{ //Insertions of non N
			//BaseRef es caracter o entero
			/*if (BaseRead[ii]==0x00) mask=0x00;
				else mask=BitsBase(BaseRead[ii], 0x00); //0,1,2,3*/
			mask=BaseRead[ii];
		}

		aux=aux|mask;
	}

	//if ((strand=='R')&&(lendesc==5)) printf("BASES %u AUXDEF %u \n", mask, aux);

	*aux_i=i;
	return(aux);
};

//Evaluates the Bits of the operation in the Forward case ()
uint8_t BitsOperF(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, int *ii ){ //puede que el strand no se necesite
	uint8_t aux=0x0, NDelC=1;
	int i=*ii;
	//OPER Y CHAR O SON ENTEROS ? pero se espera que su representacion numerica sea equivalente
	switch (oper[i]){
		case 'd': //ACTUALIZAR PARA COMPARAR CONTRA VALOR ENTERO DIRECTO
			//Calcular deleciones contiguas (4 casos), la operacion siguiente es una d y el offset es cero., hasta 4 veces

			//ojo con diferencias del STRAND
			while (((i-1)>=0)&&(offset[i]==0)&&(oper[i-1]=='d')&&(NDelC<=4)){
				NDelC++; 				i--;
			}
			if (NDelC==5) {NDelC--;i++;}
			/*while ((i<lendesc-1)&&(offset[i+1]==0)&&(oper[i+1]=='d')&&(NDelC<=4)){ //adaptation from te od one
				NDelC++; 				i++;
			}*/ //QUIT 1403014

			switch(NDelC){
				case 1:aux=0x2;break;//aux=CreateMask8B(2,1);  // Delecion Simple 0x010;
				case 2:aux=0x3;break;//aux=CreateMask8B(2,2); //0x011
				case 3:aux=0x4;break;//CreateMask8B(3,1); //0x100
				case 4:aux=0x5;break;//CreateMask8B(3,1)|CreateMask8B(1,1);  //Delecion Cuadruple Consecutiva //aux=0x101
			}
		break;
		case 'i':

			//ojo con diferencias del STRAND
			if((0x4==baseRead[i])||('n'==baseRead[i])){//Calcular inserciones en N
				//Si la operacion siguiente es una i y el offset es cero y la base del READ es la misma
				aux=0x7;
			}else{ //Insercion Simple
				aux=0x1;
			}
			break;
		case 's':
		//ojo con diferencias del STRAND
			if(((i-1)>=0)&&(offset[i]==0)&&(oper[i-1]=='s')&&(baseRead[i-1]==baseRead[i])){//Sustituciones contiguas, si operacion siguiente es s y el offset es cero y la base del READ es la misma
				//lendesc empieza en uno o a cero ?
/*OJO*/				aux=0x6;//CreateMask8B(3,2);//=0x110
				//SALTAR el numero de la operacion revisada//i=i+1;
			i--;
			}else{ //Sustitucion Simple
				aux=0x0;//Innecesario
			}
		break;
	}//Endwitch
	*ii=i;
	return (aux);

};//EndBitsOper


uint8_t BitsOperR(uint8_t *oper, uint8_t *baseRead, uint16_t *offset, uint16_t lendesc , int *ii ){ //puede que el strand no se necesite
	uint8_t aux=0x0, NDelC=1, i=*ii;
	//OPER Y CHAR O SON ENTEROS ? pero se espera que su representacion numerica sea equivalente
	switch (oper[i]){
		case 'd': //ACTUALIZAR PARA COMPARAR CONTRA VALOR ENTERO DIRECTO
			//Calcular deleciones contiguas (4 casos), la operacion siguiente es una d y el offset es cero., hasta 4 veces

			//ojo con diferencias del STRAND
			while ((i<lendesc-1)&&(offset[i+1]==0)&&(oper[i+1]=='d')&&(NDelC<=4)){
				NDelC++; 				i++;
			}
			if (NDelC==5) {NDelC--;i--;}
			switch(NDelC){
				case 1:aux=0x2; break; // Delecion Simple 0x010; CreateMask8B(2,1)
				case 2:aux=0x3; break;//0x011 CreateMask8B(2,2)
				case 3:aux=0x4; break;//0x100 CreateMask8B(3,1)
				case 4:aux=0x5; break; //Delecion Cuadruple Consecutiva //101 CreateMask8B(3,1)|CreateMask8B(1,1)
			}
		break;
		case 'i':
			//printf ("*******oper[i] %c baseRead[i] %c ",oper[i], baseRead[i]);//OJO NUMERO NO LETRA
			//ojo con diferencias del STRAND
			if((0x4==baseRead[i])/*||(110==baseRead[i])*/){//Calcular inserciones en N
				aux=0x7;
				//i++; ya no por ser ins simple
			}else{ //Insercion Simple
				aux=0x1 ; //0x001;

			}
		break;
		case 's':
		//ojo con diferencias del STRAND
			if((i<lendesc-1)&&(offset[i+1]==0)&&(oper[i+1]=='s')&&(baseRead[i]==baseRead[i+1])){//Sustituciones contiguas, si operacion siguiente es s y el offset es cero y la base del READ es la misma
				//lendesc empieza en uno o a cero ?
/*OJO*/				aux=0x6;//110;CreateMask8B(3,2)
				i++;
				//SALTAR el numero de la operacion revisada//i=i+1;
			}else{ //Sustitucion Simple
				aux=0x0;//Innecesario
			}
		break;
	}//Endswitch
	*ii=i;
	return (aux);
};//EndBitsOperR

//Distancia De la BASE en la REF a la Base en el READ, devuelve un numero del 0 al 3, al cual eventualmente tocara sumarle 1
uint8_t BitsBase(uint8_t BRead, uint8_t BRef){
	//calcula la distancia entre la base de la referencia y la base del Read
	//VectorCircular 0:A 1:C 2:G 3:T  4:N
	uint8_t aux=0x0;
	int auxInt;
	//if ((BRead!='_')&&(BRef!='_')) auxInt=abs((BaseIndex(BRead)-BaseIndex(BRef)))%5;
	if ((BRead!=0x9)&&(BRef!=0x9)){ //MUST BE Differente from dash '-'

		if (BRead>BRef){
			//auxInt=abs( ((int)BaseIndex(BRead)-(int)BaseIndex(BRef)) );
			auxInt=abs((int)BRead-(int)BRef);
		}else{
			if (BRead==BRef) printf("Error Grave entre bases iguales Base1 %u Base2 %u \n", BRead,  BRef);
			else{
				auxInt=((5-(int)BRef)+(int)BRead);//%5;//No sumar el 1 sumarlo al calclular Index2Base
				//auxInt=((5-(int)BaseIndex(BRef))+(int)BaseIndex(BRead));//%5;//No sumar el 1 sumarlo al calclular Index2Base
				//aux=((int)BaseIndex(BRead)+(int)BaseIndex(BRef)+1)%5;//No sumar el 1 sumarlo al calclular Index2Base
			}
		}
		//printf("Diferencias Ref %i   Read %i  Distancia %i \n", BRef,BRead, aux+1);
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
	}//endif
	//printf("Diferencias Ref %i  Read %i  Distancia %i auxInt %i\n", BRef,BRead, aux+1, auxInt);
	return(aux);
};

uint8_t BaseIndex(char Base){   	//CAlculo del indice de la Base
	switch(Base){
		case 'A':return(0);break;
		case 'a':return(0);break;
		case 'C':return(1);break;
		case 'c':return(1);break;
		case 'G':return(2);break;
		case 'g':return(2);break;
		case 'T':return(3);break;
		case 't':return(3);break;
		case 'N':return(4);break;
		case 'n':return(4);break;
		case '_': printf("___"); return(8);
		break;
		default: printf("Error Grave de Base NO RECONOCIDA %c \n",  Base );
	}
}

char Index2Base(uint8_t Index){   	//CAlculo del indice de la Base
	switch(Index){
		case 0:return('A');break;
		case 1:return('C');break;
		case 2:return('G');break;
		case 3:return('T');break;
		case 4: case 78: return('N');break;
		default: printf("Index unknown %i \n",  Index );
	}
}

/*8 bits menos significativos del offset. Los 2 bits + significativos en el 3er Entero, en las posiciones menos significativas*/
uint8_t Offset( uint16_t offset, uint8_t *rest){
	uint8_t mask=0x01, aux=0x0, auxRest=0x0, mod, bitCnt=0;

    while (offset>0){

    	if (bitCnt<8){
    		mod=offset%2; offset=offset/2;
	    	if (mod==1) aux=aux|mask;
	    	mask=mask<<1;  //HAcia la izquierda porque la modificacion ocurre en la mascara y no el aux
	    	bitCnt++;
	    	if (bitCnt==8) mask=0x01;
    	}else{
    		mod=offset%2; offset=offset/2;
    		if (mod==1) auxRest=auxRest|mask;
    		mask=mask<<1;
    	}
    }//endwhile

    if (offset==1) //Digito del ultimo cociente
    	if (bitCnt>=8){
    		auxRest=auxRest|mask;
    	}else{
    		aux=aux|mask;
    	}

    //construir por separado el resto para que se pueda pocesar en la siguiente instruccion
    *rest=auxRest;
	return(aux);
};

//¿Como resolver asunto del more frags?, traer esa info del alineador ?
uint8_t Preambulo(uint8_t moreFrags, char strand, uint16_t  lendesc, uint16_t EdDis ){
	uint8_t mask=0x01, aux=0x0;

	//1. Considerar caso de dos preambulos perfectos contiguos. si hay espacio?,
		//Como si solo tengo la info de un read al momento?, pero lo mismo ocurre con el mostfrags!
	if (moreFrags==1) { //MoreFrags *RESOLVER
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
		if (strand=='F') mask=0x4;//CreateMask8B(3,1);//0x100;   //Forward
		if (strand=='R') mask=0x5;//CreateMask8B(3,1)|CreateMask8B(1,1);//0x101;   //Reverse
		if (strand=='C') mask=0x6;//CreateMask8B(3,2);//0x110;   //Complement
		if (strand=='E') mask=0x7; //CreateMask8B(3,3);//0x111;   //rEvErsE complEmEnt
		//if (strand=='T') mask=0x011; //another Transform
	};
    if ((lendesc==1)&&(EdDis==1)){
			if (strand=='R') mask=0x5;//CreateMask8B(3,1)|CreateMask8B(1,1);//0x101;   //Reverse
			if (strand=='E') mask=0x7; //CreateMask8B(3,3);//0x111;   //rEvErsE complEmEnt
    }
	aux=mask|aux;
	return(aux);
};
/*******************************************************************************************************************************/
/************************************************END of Compression ( Instr to Binary  )*****************************************/
/*******************************************************************************************************************************/



/********************Reconstruccion de la Instruccion (Previamente probadas) **************/
  /*Este modulo pretende reconstruir y verificar el read en funcion de las instrucciones generadas*/
  int VerificRead(Olg olg, Instr *MInst, char *ref1, int lenR, uint8_t Match, uint16_t *Of, unsigned char *Op, uint8_t *Ba){
  	int tam,w,iguales=0, posRead, posRef, ofAct,band=1, borradasF=0, borradasR=0, InsR=0, InsF=0;
  	static int raros=0;
  	char *auxR, *MiRead;
  	FILE *aF, *aU;

  	//1.Calcular la longitud del read a generar y guardarla en tam:

  	if (olg.strand=='R'){
  		if (olg.NDel>0){
  			tam=olg.lonRead+olg.NDel;
  		}else{
  			tam=olg.lonRead;
  		}
  	}else{
  		if (olg.NDel>0){
  			tam=olg.lonRead+olg.NDel;
  		}else{
  			tam=olg.lonRead;
  		}
  	}

  	//CAlculo mi Read Inicial a partir de la referencia
  	auxR=ref1+olg.loc;
  	if (MiRead=(char*)malloc(sizeof(char)*((tam)+5))){
  	 	MiRead[0]='\0';
  		strncat(MiRead, auxR, tam);
  	}else
  		printf ("No hay memoria suficiente para crear MiRead");

  	//Solo usado para verificacion de ReadsRaros
  	if (aF=fopen("Verifica.txt","a")){
  		fprintf (aF,"\n /*****************************PRINT INI PAR READ %i **************************************/", olg.NREAD);
  			fprintf (aF,"\n %s \n",MiRead);fflush(stdout);
  		if (aF) fclose(aF);
  	}

      //Aplicar Operaciones por cada instruccion.
  	//Regla : las ausencias de enteros en el desc se consideran un cero.
  	if (Match=='F') { //F
  				//Si es FORWARD:  evalúa el descriptor de izquierda a derecha, se ignora el ultimo entero del descriptor (el de mas a la izq) .  Se parean los números con el carácter de la izquierda. La primera instrucción en la ED se ejecuta de primera de derecha a izq sobre la REF.

  		posRead=0;
  		for (int y=olg.lendesc-1;y>=0 ;y--){ //cada operacion, pendiente con edit distance
  			char aux=DesconvertBase(Ba[y],1);

  			switch(Op[y]) {
  			   case 's': // De acuerdo al MD sustitucion en en read, -->sustitucion en la referencia

  				   posRead=posRead+Of[y+1]; //Calculo la posicion del Offset en el read
  				   //printf ( "         PosRead % i       ", posRead);
  				   //printf (aF,"\n Instruccion %i dice meta %c letra %c en la pos %i offset %i\n", y,Op[y],aux,posRead, Of[y]);
  				   if ((aux!='2')&&(aux!='3')&&(aux!='4')){
  					   MiRead[posRead]=aux;
  					   posRead++;			  //Descuento la base actualizada
  				   }else{   //Caso de Varias N seguidas

  					  /*int tope=aux - '0';  //Dice cuantas N seran sutituidas
  					  // printf ("\n Aux %c tope %i \n", aux, tope);
  					   for (int b=0; b<tope;b++){
  						   MiRead[posRead]='N';   //este tambien podia ser, de hecho deberia ser porque las n seguidas van en funcion del read.
  						   //MiRead[posRead]=olg.read[posRead];//aux=ref1+olg.loc //En el caso real debe ser por NN
  	   					   posRead++;			  //Descuento la base actualizada
  					   }*/
  				   }
  				   //printf ("\n Modificada la posicion %i con la letra %c y el vector %c \n", posRead, aux, MiRead[posRead+1+Of[y]]);

  			   break;
  			   case 'd': // De acuerdo al MD delecion en la referencia -->insercion en el read al final
  			  		posRead=posRead+Of[y+1]; //Calculo la posicion del Offset en el read

  					if ((aux!='2')&&(aux!='3')&&(aux!='4')){ //Nucleotidos simples
  					   for(int q=posRead;q<tam-1;q++){
  						  MiRead[q]=MiRead[q+1];
  					    }
  					    borradasF++;
  						//no se decrementa ya que la pos igue siendo en POSREAD
  					}

  			   break;
  			   case 'i': // De acuerdo al MD una insercion -->delecion
  				   //Inserte en el read para que se asemeje a la referncia, corresponde a eliminar en la referencia (MiRead)
  				   //la insercion ocurre en l a POS SIGUIENTE, no en la exacta

  				   //printf ("\n Instruccion %i dice Oper %c letra %c en la posRead %i offset %i\n", y,Op[y],aux,posRead, Of[y]);
  				   posRead=posRead+Of[y+1]; //Calculo la posicion del Offset en el read

  				   if ((aux!='2')&&(aux!='3')&&(aux!='4')){ //Nucleotidos simples
  					   for(int q=(tam-1);q>posRead;q--){
  						  MiRead[q]=MiRead[q-1];
  						}

  						MiRead[posRead]=aux;
  						posRead++;
  				   }
  						/*else{// CONSIDERA inserciones con varias n.
  					   int tope=aux - '0';  //Dice cuantas N seran insertadas
  					   //printf ("\n Aux %c tope %i \n", aux, tope);
  					   for(int q=(olg.lonRead-1);q>=posRead+tope+1;q--) { //Igualmente cuidar el sentido del desplaczamiento
  						   MiRead[q]=MiRead[q-tope];
  					   }
          			   for (int b=0; b<tope;b++){
          				   MiRead[posRead+b+1]='N';
  						}
          			   posRead++;
    					 //}//endelse

  					 //printf ("\n Modificada la posicion %i con la letra %c y el vector %c \n", posRead, aux, MiRead[posRead+1+Of[y]]);
  					*/
  				   break;
  				   case '_': case '?': case 'X'://posRead=posRead-Of[y];//Solo informacion de offset
  				   break;
  			}//end switchF

  		}//endForyF

  	}else{  	//REVERSE Match==1
  		//Si es Reverse: Se evalúa el descriptor de derecha a izquierda, se ignora el primer entero del descriptor (el de mas a la izq) .  Se parean los números con el carácter de la izquierda. La primera instrucción en la ED se ejecuta de primera de derecha a izq sobre la REF.
          posRead=(tam-1)-olg.NIns;

  		for (int y=0;y<olg.lendesc;y++){

  			//Solo para correccion de errores
  			/*if (aU=fopen("PosRead.txt","a")){
  			fprintf (aU,"\n READ %i PosREAD %i instruccion %i tam %i \n", olg.NREAD,posRead, y, tam);fflush(stdout);
  			fclose(aU);
  			}*/

  			char aux=DesconvertBase(Ba[y],1);

  			switch(Op[y]) {
  			   case 's': // De acuerdo al MD sustitucion en read, -->sustitucion en la referencia

  				   posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read
  				   /*if (aF=fopen("Verifica.txt","a")){
  				   fprintf (aF,"\n Instruccion %i dice meta %c letra %c en la pos %i offset %i\n", y,Op[y],aux,posRead, Of[y]);
  				   fclose(aF);
  				   }*/

  				   if ((aux!='2')&&(aux!='3')&&(aux!='4')){
  					   MiRead[posRead]=aux;
  					   posRead--;			  //Descuento la base actualizada
  				   }
  				   /*else{ //sustitucion de Varias N seguidas
  					   int tope=aux - '0';  //Dice cuantas N seran sutituidas
  					   printf ("\n Aux %c tope %i \n", aux, tope);
  					   for (int b=0; b<tope;b++){
  						   //MiRead[posRead]='N';   este tambien podia ser, de hecho deberia ser porque las n seguidas van en funcion del read.
  						   MiRead[posRead]=olg.read[posRead];//aux=ref1+olg.loc //En el caso real deberia ser por N
  	   					   posRead--;			  //Descuento la base actualizada
  					   }
  				   }*///EndElse

  			   break;
  			   case 'd': // De acuerdo al MD delecion en la referencia -->insercion en el read
  				   posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read
  				   /*if (aF=fopen("Verifica.txt","a")){
  				   	fprintf (aF,"\n Instruccion %i dice meta %c letra %c en la pos %i offset %i\n", y,Op[y],aux,posRead, Of[y]);
  				   	fclose(aF);
  				   }*/
  				  if ((aux!='2')&&(aux!='3')&&(aux!='4')){ //Nucleotidos simples
    				     for(int q=posRead;q<tam-1;q++){ /***********ver si este tope con strlen o incluir NINS************/
  					    MiRead[q]=MiRead[q+1];
  				     }
  				  }
  				  posRead--;
  				  borradasR++;
  			   break;

  			   case 'i': // De acuerdo al MD una insercion -->delecion
  				   //Inserte en el read para que se asemeje a la referncia, corresponde a eliminar en la referencia (MiRead)
  				   //la insercion ocurre en l a POS SIGUIENTE, no en la exacta

  				   //printf ("\n Instruccion %i dice Oper %c letra %c en la posRead %i offset %i\n", y,Op[y],aux,posRead, Of[y]);
  				   posRead=posRead-Of[y]; //Calculo la posicion del Offset en el read
  				   /*if (aF=fopen("Verifica.txt","a")){
  				   fprintf (aF,"\n Instruccion %i dice meta %c letra %c en la pos %i offset %i\n", y,Op[y],aux,posRead, Of[y]);
  				   fclose(aF);
  				   }*/

  				   if ((aux!='2')&&(aux!='3')&&(aux!='4')){ //Nucleotidos simples
  					   for(int q=(tam-1);q>posRead+1;q--){   //utilizar la longitud de la variable mas quela del read
  						  MiRead[q]=MiRead[q-1];
  					   }
  					   MiRead[posRead+1]=aux;  //.Puede que toque sumar 1

  				   } /*else{// CONSIDERA QUE HAY inserciones con varias n.
  					   int tope=aux - '0';  //Dice cuantas N seran insertadas
  					   //printf ("\n Aux %c tope %i \n", aux, tope);

  					   for(int q=(olg.lonRead-1);q>=posRead+tope+1;q--) {
  						   MiRead[q]=MiRead[q-tope];
  					   }
          			   for (int b=0; b<tope;b++){
          				    MiRead[posRead+b+1]='N';
  					   }
          			   	 //posRead=posRead--; //ojo
  				     }*///end else

  				   InsR++;
  				   //printf ("\n Modificada la posicion %i con la letra %c y el vector %c \n", posRead, aux, MiRead[posRead+1+Of[y]]);

  				   break;

  			   case '_': case '?': case 'X': case 'x'://posRead=posRead-Of[y];//Solo informacion de offset
  			   break;

  			}//end switch

  		}//endFory

  	}//else Reverse
  	MiRead[olg.lonRead]='\0';

  	//Verificacion de la reconstruccion
  	if ((band)&&(strcmp(olg.read,MiRead)!=0)){
  		if (aF=fopen("Verifica.txt","a")){
  			fprintf (aF,"/Original***************************PRINT End PAR*********************************MiRead/ \n");
  			fprintf (aF,"\n Caso raro de error en el read numero : %i \n",olg.NREAD);

  			fprintf (aF,"\n %s\n",olg.read);
  			fprintf (aF,"\n %s\n",MiRead);
  			raros++;
  			printf("\n TOTAL RAROS %i, read %i, sentido %c, NiNs %i, Ndels %i y posRead %i editD %i ",  raros ,olg.NREAD, olg.strand,  olg.NIns, olg.NDel, posRead, olg.EditDist);
  			if (aF) fclose(aF);
  		}
  	}
  	return (strcmp(olg.read,MiRead));//Compara el READ generado con el del olg
  };

  //Desconvierte un numero (codigo) a una base en caracter
  char DesconvertBase(uint8_t Ba, int repeat){
  	char sal;
  	switch(Ba) {
  	   case 0x0: sal='A';
  	   break;
  	   case 0x1: sal='C';
  	   break;
  	   case 0x2: sal='G';
  	   break;
  	   case 0x3: sal='T';
  	   break;
  	   case 0x4: sal='N';
  	   break;
  	   /*case 5:   sal='2';//NN
  	   break;
  	   case 6:   sal='3';//NNN
  	   break;
  	   case 7:   sal='4'; //NNNN
  	   break;*/
  	   default : sal='_'; //code for Error
  	 }
  	return (sal);
  };

  //Convierte las bases involucradas en enteros para meterlos en las Estructuas dinamicas
  uint8_t ConvertBase(int repeat, char rd_ltr){	//Corresponden a la OPER 0=A, 1=C, 2=G, 3=T, 4=N, 5=NN, 6=NNN, 7=sNNNN
  uint8_t sal;
  rd_ltr=toupper(rd_ltr);
 // if (repeat<=1){
  	switch(rd_ltr) {
  	   case 'A': case 'a': sal=0x0;
  	   break;
  	   case 'C': case 'c': sal=0x1;
  	   break;
  	   case 'G': case 'g': sal=0x2;
  	   break;
  	   case 'T': case 't': sal=0x3;
  	   break;
  	   case 'N': case 'n': sal=0x4;
  	   break;
  	   default : sal=0x9; //printf ( "************ESTA BASE DA ERROR ********** %c", rd_ltr );
  	}
  	//if (repeat == 0) sal=100;
/*  }else{
  	switch(repeat) {
  		case 2:   sal=5;//NN;
  		break;
  		case 3:   sal=6;//NNN
  		break;
  		case 4:   sal=7; //NNNN
  		break;
  		default : sal=100; //code for Error
  	 }
  }//else;
*/
  return (sal);
  };

  //it separates the read, the location, the strand and the match descriptor from the input BAM object
  Olg MIparse_input(bam1_t *b, char *ref, int ref_len, long *contU){
  	Olg ret_var;
  	int i,TamRead;

  	bam1_core_t core = b->core;
  	uint8_t *seq = bam1_seq(b);

  	static char *read_s=NULL;
  	static char *qual_s=NULL;
  	int max_bytes=10000;

  	if (read_s==NULL)  		        read_s=(char*)malloc(sizeof(char)*max_bytes);

  	if(core.l_qseq >= max_bytes)  	printf ("Unsupported l_gseq, exceed the allowed amount of bytes");

  	for (i=0; i<core.l_qseq; i++) { ///***********************************LECTURA DEL READ ********************************///

  		int base = bam1_seqi(seq, i);  //to get the base

  		switch (base) {                //To assign each base
  		case 1:
  			read_s[i] = 'A';
  			break;
  		case 2:
  			read_s[i] = 'C';
  			break;
  		case 4:
  			read_s[i] = 'G';
  			break;
  		case 8:
  			read_s[i] = 'T';
  			break;
  		case 15:
  			read_s[i] = 'N';
  			break;
  		default:
  			//TODO: should never get here...throw exception?
  			read_s[i] = 'X';
  			break;
  		}
  	}
  	read_s[core.l_qseq] = '\0';  // Null terminate my own strings

  	char strand;
  	if ((core.flag & 0x10) == 0x10)
  		strand = 'R';
  	else
  		strand = 'F';  //cuidado porque no es el unico caso dentro del SAM, aunque sin duda parecieran los unicos importantes

  	if ((core.flag & 0x4) == 0x4) {*contU=*contU+1; strand='U';} //Conteo de unmmaped
  	int Edist=0;

  	/******************************Si quiero la data en objetos***************************/
  	ret_var.read = read_s;
  	ret_var.qvals = qual_s;
  	ret_var.loc = core.pos+1;
  	ret_var.strand = strand;
  	int Ni, Nd;
  	ret_var.descriptor = create_match_descriptor(b, ref, ref_len, read_s, strand, &TamRead, &Edist, &Ni, &Nd);

  	ret_var.NDel=Nd;
  	ret_var.NIns=Ni;
  	ret_var.lonRead=TamRead;
      ret_var.EditDist=Edist;
      return ret_var;

  }

  void add_new_instruction(Instr *instr_ar, int i, char type, char ref_ltr, char rd_ltr,int repeat, int offset){
  	if(i>=INSTR_LENGTH){
  		fprintf(stderr,"FATAL ERROR!!!Instr index greater than total supported instructions\n");
  		exit(0);
  	}

  	/******************************Si quiero la data en objetos***************************/
  	instr_ar[i].offset=offset;
  	instr_ar[i].type=type;
  	instr_ar[i].ref_ltr=ref_ltr;
  	instr_ar[i].rd_ltr=rd_ltr;
  	instr_ar[i].repeat=repeat;

  	/***********************Si quiero la data en Arreglos Dinamicos***********************/
  	/*colocar en variables las constantes almacenadas.*/
  }


  /* In: Info de alineacion. //Data of *len will contain the length of that array.
   *
   * El objetivo del estructurador de instrucciones radica en procesar los resultados provenientes de la alineacion
   * Calcula el prefijo de cada Instrucción, modifica el mapa, y genera el arreglo de mapeo considerando los offsets, y OPCODE.
   * Los deja concatenados. Este trabaja sobre una solo READ a la vez.
   *
   * En principio se calculan los offsets de todas las instrucciones y se determina
   * la cantidad de bits necesarias para cada instruccion.
   * offst_bits es el numero de bts que reresenta un  offset , rest_bits es el resto de bits de la isntruccion.
   *
   * Out: Retorna el array de instrucciones que sicroniza la informacion del Cigar con la del MD
   * */

  Instr *parse_line(Olg olg, int *len, char *ref){ //Tamaño del offset en cifras decimales (entero corto o arreglo de binario?) Deben ser 10 bits, para 1088 posiciones
  	int offst[100];      			//Cantidad de offsets para este read
  	char *dscr=olg.descriptor;		//Lista de descriptor de las operaciones obtenidas poe el alineador
  	char *rd=olg.read;				//En el caso real, yo tendre los reads porque vienen del FASTQ.
  	Instr *ret_ar=(Instr*)malloc(INSTR_LENGTH*sizeof(Instr));    //Donde voy ubicando la instruccion a guardar
  	int abs_offst_rd=strlen(rd)-1,abs_ref_rd=strlen(rd)-1 ;  //Lo inicio en la longitud del Read-1
  	int instr_len=0;              	//Para ir contando la longitud de la isntruccion
  	int offst_len=0, num_len=4;		//Para contar la longitud de este offset,
  	int i=0, longRCalc=0;			//Longitud de Read calculada en base a los offsets
      //int  Ndel=0, Nin=0;

  	char prev='0'; //'A' if nucleotide, 'n' if numeric, '0' if at the beginning, 'd' for indel.
  	char num[5]={'0','0','0','0','\0'};    //Donde almaceno el offset actual
  	int last_offst=0;					   //Indica la posicion del offset mas reciente
  	int open_indel=0, total_ins=0, cnt=0;  //Indica si hay una operac de indel abierta, y cuenta cuantas ins, y general

  	for (i=strlen(dscr)-1;i>=0;i--){ 	   //Se evalua el descrptor de derecha a izquierda

  		if (dscr[i]>='0' && dscr[i]<='9'){ // el elemento procesado es numerico

  			if (prev!='n'){       	//Si el elemento anterior NO es numerico (este es el primer numero de este grupo)
  				num_len=4;        	//Inicializo el tamaño de este numero, maximo 5 digitos para 99999
  				num[0]='0';num[1]='0';num[2]='0';num[3]='0'; num[4]='\0';   //Los inicializo
  			}
  			num[num_len--]=dscr[i]; //Asigno el numero leido en el lugar correspondiente
  			prev='n';               //para guardar que el ultimo elemento procesado es un numero

  			if (i==0){  			   //para guardar el ultimo numero que permite completar el read
  				last_offst=atoi(num);  //Actualizo el ultimo offset //convertindiendlo en entero
  				longRCalc= longRCalc+last_offst;
  				//Meterlo en una instruccin vacia y Nula
  				add_new_instruction(ret_ar, instr_len++, '_', '_', '_', 0, last_offst);
  			}
  		}
  		else if(dscr[i]>='A' && dscr[i]<='T'){ // el elemento procesado es un NUCLEOTIDO

  			if(prev=='0'){                     //Si se esta iniciando este grupo
  				offst[offst_len++]=0;          //AGrego uno mas al arreglo de offsets y lo inicializo
  				last_offst=0;                  //Inicializo el offset anterior
  				//Guardar la instruccion

  				if (abs_offst_rd>=0)  add_new_instruction(ret_ar, instr_len++, 's', dscr[i], rd[abs_offst_rd], 1, last_offst);
  				longRCalc++;
  				abs_offst_rd--;                //disminuyo 1 al valor de la posicion del offset read evaluado,
  			}else
  				if(prev=='n'){          	   //Si El ultimo elemento procesado es un NUMERO

  					if(open_indel){            //Si hay una operacion de insercion eliminacion abierta
  						fprintf(stderr, "FATAL ERROR!!!Cannot have both insertions and deletions in one indel block\n");
  						exit(1);
  					}

  					offst[offst_len++]=atoi(num)+1;  //Agrego a la lista de offsets el recientemente calculado
  					last_offst=atoi(num);//+1;       //Actualizo el ultimo offset
  					longRCalc=longRCalc+atoi(num);

  					abs_offst_rd=abs_offst_rd-(last_offst); //Actualizo la nueva posicion offset del Read

  					if (abs_offst_rd>=0) add_new_instruction(ret_ar, instr_len++, 's', dscr[i], rd[abs_offst_rd], 1, last_offst);

  					abs_offst_rd--; 	       		//Desplazo a la siguiente posicion de offset de read.
  					longRCalc++;
  				}else{  //(prev!=n)

  					if(open_indel){         		//caso de la delecion
  						if (abs_offst_rd>=0) add_new_instruction(ret_ar, instr_len++, 'd', dscr[i], rd[abs_offst_rd], 1, last_offst);
  						last_offst=0;      //Reinicio el offset previo
  						//Ndel++;
  						prev='A';		   //Marco el previo como nucleotrido para que se retome un nuevo error
  						longRCalc++;       //Por cada delecion se asume que en el read hay un elemento mas que en la referencia.
  						continue;
  					}
  					offst[offst_len++]=0;
  					last_offst=0;
  					abs_offst_rd=abs_offst_rd-last_offst;

  				/*	if ((abs_offst_rd>=0) &&(ret_ar[instr_len-1].repeat<=3)&&(rd[abs_offst_rd]=='N' && ret_ar[instr_len-1].rd_ltr=='N' && ret_ar[instr_len-1].type=='s'))
  					{ 	//Caso de susitucion Varias N consecuctivas
  						ret_ar[instr_len-1].repeat+=1; longRCalc++;//Contador de repeticiones de N
  						abs_offst_rd--;//recien agregado el 061016
  					}else{*/
  						if (abs_offst_rd>=0) add_new_instruction(ret_ar, instr_len++, 's', dscr[i], rd[abs_offst_rd], 1, last_offst);
  						abs_offst_rd--;longRCalc++;
  				//	};
  				} //ende else (prev==n)
  			prev='A';				//Marco el previo como  un nucleotido

  			if ((i==0)&&(olg.strand=='F')){     //ultimo caracter y Sentifo F, y ya sabiamos q era nucleotido
  				add_new_instruction(ret_ar, instr_len++, '_', '_', '_', 0, 0);
  				longRCalc++;
  			}

  		}else
  			if(dscr[i]=='$'){  		//Inicio del Segmento de Indel*
  				open_indel=1;
  				if(prev=='0'){      //SiEstoyalInicio
  					last_offst=0;
  				}else
  					if(prev=='n'){  //Sielanterior es un numero
  						last_offst=atoi(num);
  					    abs_offst_rd=abs_offst_rd-(last_offst);
  					}else{ 			//bp or another indel is preceding
  						last_offst=0;
  					}
  					prev='d';       //Marco el previo como un indel
  			}else
  				if(dscr[i]=='^'){ 	//Fin del Segmento*..
  					open_indel=0;
  					if(prev=='n'){
  						total_ins=atoi(num);
  						for(cnt=0;cnt<total_ins;cnt++){

  							if (abs_offst_rd>=0){

  								/*if((rd[abs_offst_rd]=='N' &&(ret_ar[instr_len-1].repeat<=3) &&  ret_ar[instr_len-1].rd_ltr=='N' && ret_ar[instr_len-1].type=='i')){ //	insercion de varias N
  									ret_ar[instr_len-1].repeat+=1; //Contador de repeticiones de N INSERCION
  									abs_offst_rd--;
  	   							    //No cuento crecimiento del read
  									//Notese que esto SI Esta programado pensando en las N del read y no de la ref
  								}else {*/
  									add_new_instruction(ret_ar, instr_len++, 'i', ref[abs_offst_rd], rd[abs_offst_rd], 1, last_offst); //dsrc[i en ref]
  									longRCalc=last_offst+longRCalc; //Cuento el offset
  								//};

  								last_offst=0;
  								abs_offst_rd--;
  								//Al insertar en la referencia, el tamaño del read no cambia
  							}
  						}//for cnt
  					}//If prev==n
  					prev='0';
  				}//IF dscr[i]=='^'

  	}//END FOR
  	(*len)=instr_len; 			//Cantidad de elementos en el descriptor

  	return ret_ar;
  }

  /***************************************************QSORT********************************************/
  void swap (uint32_t *arr, uint32_t a, uint32_t b){
  	uint32_t aux;
  	aux=arr[b]; arr[b]=arr[a]; arr[a]=aux;
  }

void insertionSortInd( uint32_t *arr, int low, int high, uint32_t *Indexes)  {
  int j;
  uint32_t key, aux;
  for (int i = low + 1; i <= high; i++){
    key = arr[i];
    aux=Indexes[i];
    j = i-1;

    /* Move elements of arr[low..i-1], that are greater than key, to one position ahead  of their current position */
    while (j >= low && arr[j] > key){
      Indexes[j+1]=Indexes[j];
      arr[j+1] = arr[j];
      j = j-1;
    }
    Indexes[j+1]=aux;
    arr[j+1] = key;
  }
}

  /* Function to sort an array using insertion sort*/
  void insertionSort(uint32_t *arr, int low, int high)
  {
     int key, j;
     for (int i = low + 1; i <= high; i++)
     {
         key = arr[i];
         j = i-1;

         /* Move elements of arr[low..i-1], that are
            greater than key, to one position ahead
            of their current position */
         while (j >= low && arr[j] > key)
         {
             arr[j+1] = arr[j];
             j = j-1;
         }
         arr[j+1] = key;
     }
  }


int randomizedPartitionInd ( uint32_t *arr, int low, int high, uint32_t *Indexes)	  {
// find random element
  int rand_pivot = low + rand() % (high - low + 1);
  // swap random pivot to end.
  swap(arr, rand_pivot, high);
  swap(Indexes, rand_pivot, high);

  int pivot = arr[high];    // pivot
  int i = (low - 1);  // Index of smaller element

  for (int j = low; j <= high- 1; j++)   {
  // If current element is smaller than or   // equal to pivot
	  if (arr[j] <= pivot)      {
		  i++;    // increment index of smaller element
		  swap(arr,i, j);
		  swap(Indexes, i, j);
      }
   }
   swap(arr,i + 1,high);swap(Indexes,i + 1,high);
   return (i + 1);
}




  /* This function takes random element as pivot, places
     the pivot element at its correct position in sorted
     array, and places all smaller (smaller than pivot)
     to left of pivot and all greater elements to right
     of pivot */
  int randomizedPartition (uint32_t *arr, int low, int high)
  {
  	// find random element
  	int rand_pivot = low + rand() % (high - low + 1);
  	// swap random pivot to end.
      swap(arr, rand_pivot, high);

      int pivot = arr[high];    // pivot
      int i = (low - 1);  // Index of smaller element

      for (int j = low; j <= high- 1; j++)
      {
          // If current element is smaller than or
          // equal to pivot
          if (arr[j] <= pivot)
          {
              i++;    // increment index of smaller element
              swap(arr,i, j);
          }
      }
      swap(arr,i + 1,high);
      return (i + 1);
  }



  void optimizedQuickSortIndex( uint32_t *Pos, int low, int high, uint32_t *Indexes)
  {
    	while (low < high)
    	{
    		// do insertion sort if 10 or smaller
    		if(high - low < 10)
    		{
    			insertionSortInd(Pos, low, high, Indexes);
    			break;
    		}
    		else
    		{
    			int pivot = randomizedPartitionInd(Pos, low, high, Indexes);
    			// perform tail call optimizations
    			// recurse on smaller sub-array
    			if (pivot - low < high - pivot) {
    				optimizedQuickSortIndex(Pos, low, pivot - 1,Indexes);
    				low = pivot + 1;
    			} else {
    				optimizedQuickSortIndex(Pos, pivot + 1, high,Indexes);
    				high = pivot - 1;
    			}
    		}
    	}
    }


  /* The Non-Optimized QuickSort
    arr[] --> Array to be sorted,
    low  --> Starting index,
    high  --> Ending index */
  void NonOptimizedQuickSort(uint32_t arr[], int low, int high)
  {
      if (low < high)
      {
          /* pi is randomizedPartitioning index, arr[p] is now
             at right place */
          int pi = randomizedPartition(arr, low, high);

          // Separately sort elements before
          // randomizedPartition and after randomizedPartition
          NonOptimizedQuickSort(arr, low, pi - 1);
          NonOptimizedQuickSort(arr, pi + 1, high);
      }
  }

  /* The Optimized QuickSort
    arr[] --> Array to be sorted,
    low  --> Starting index,
    high  --> Ending index */
  void optimizedQuickSort(uint32_t *A, int low, int high)
  {
  	while (low < high)
  	{
  		// do insertion sort if 10 or smaller
  		if(high - low < 10)
  		{
  			insertionSort(A, low, high);
  			break;
  		}
  		else
  		{
  			int pivot = randomizedPartition(A, low, high);
  			// perform tail call optimizations
  			// recurse on smaller sub-array
  			if (pivot - low < high - pivot) {
  				optimizedQuickSort(A, low, pivot - 1);
  				low = pivot + 1;
  			} else {
  				optimizedQuickSort(A, pivot + 1, high);
  				high = pivot - 1;
  			}
  		}
  	}
  }

  /*************************Funcion para optimizar la representacion del mapa general de correspondencias******************/
  //Toma un mapa optimizado y devuelve la representacion vectorial del mapa a 1 bit por posicion
//Previously generate data structures
  //Tome un vector de LenOpMap de 64 bits y coloque cada bit en un vector

  	//TOma un mapa que tiene lREf Cantidad de Caracteres, y devuelve un vector de enteros de 64bits con el nuevo mapa
  //Creates a binary string which represents the whole map , ready to be compressed (low level) and then stored
  void DeConstructOptimizedMap(uint8_t *MAPA , uint64_t  *MAPAO, int lRef, int LenOpMap){
  	uint64_t mask=0x01, aux, posCount=1, aux2;

  	mask=mask<<63;
  	/*DEconstructing el mapa optimizado*/
  	for  (long y=0;y<LenOpMap;y++){
  		aux=MAPAO[y];

  		 for (uint8_t bitCount=0;bitCount<64;bitCount++){
  			 //Mascara de 64 bits
  			 aux2=aux&mask;
  			 MAPA[posCount]=aux2>>63;
  			 aux=aux<<1;
  			 posCount++;
  		 }
  	}//For long y
  	//printf ();
  	//No guardo la longitud del mapa porqu elo controo con la lenRef
} /*Deconstruct End*/


  /*************************Funcion para optimizar la representacion del mapa general de correspondencias 8BITS******************/
   //Toma un mapa optimizado y devuelve la representacion vectorial del mapa a 1 bit por posicion
 //Previously generate data structures
   //Tome un vector de LenOpMap de 64 bits y coloque cada bit en un vector

   	//TOma un mapa que tiene lREf Cantidad de Caracteres, y devuelve un vector de enteros de 64bits con el nuevo mapa
   //Creates a binary string which represents the whole map , ready to be compressed (low level) and then stored
   void DeConstructOptimizedMap8(uint8_t *MAPA , uint8_t  *MAPAO, int lRef, int LenOpMap){
   	uint64_t mask=0x01, aux, posCount=1, aux2;
   	int bits =8;

   	mask=mask<<(bits-1);
   	/*DEconstructing el mapa optimizado*/
   	for  (long y=0;y<LenOpMap;y++){
   		aux=MAPAO[y];

   		 for (uint8_t bitCount=0;bitCount<bits;bitCount++){
   			 //Mascara de 64 bits
   			 aux2=aux&mask;
   			 MAPA[posCount]=aux2>>(bits-1);
   			 aux=aux<<1;
   			 posCount++;
   		 }
   	}//For long y
   	//printf ();
   	//No guardo la longitud del mapa porqu elo controo con la lenRef
 } /*Deconstruct End*/


  /*************************Funcion para optimizar la representacion del mapa general de correspondencias******************/
  	//TOma un mapa que tiene lREf Cantidad de Caracteres, y devuelve un vector de enteros de 64bits con el nuevo mapa
  //Creates a binary string which represents the whole map , ready to be compressed (low level) and then stored
  //starts at position ZERO
  void optimizeMap(uint8_t *MAPA , uint64_t  *MAPAO, int lRef, int *LenOpMap){
  	unsigned int posRel=1, posVec=0;
  	uint64_t mask=0x01, aux=0x0;
  	/*Generando el mapa optimizado*/
  	for  (int y=1;y<lRef;y++){  //ojo diferencias reerse y forwad

  		if (aux!=0) aux=aux<<1;
  		if(MAPA[y]==1){ //Meta un 1
  			aux=mask|aux;

  		}

  		if ((posRel>=64)||(y>=lRef)) {
  			if (y>=lRef) aux=aux<<(64-posRel); //o -1?
  			MAPAO[posVec]=aux;
  			//if (y<=65) printf ("  iter  %i AUX %"PRIu64 "\n",y, aux);
  			//if (y<=65) printf ("  iter  %i AUX %lu\n",y, aux);//both are the same

  			if (y<lRef){
  				 posRel=1; 			posVec++;   			aux=0x0;

  			}
  		} else {
  			posRel++;
  		}

  	}//For Int y

  	*LenOpMap=posVec;
  } /*Fin del mapa optimizado*/


  /*************************Funcion para optimizar la representacion del mapa general de correspondencias 8 bits ******************/
    	//TOma un mapa que tiene lREf Cantidad de Caracteres, y devuelve un vector de enteros de 64bits con el nuevo mapa
    //Creates a binary string which represents the whole map , ready to be compressed (low level) and then stored
    //starts at position ZERO.
  //Esta version de 8 bits pretende adaptarse a requerimientos de pbzip2 y eliminar asi overhead
  //LenOpMap8 starts in 0, so the REAL SIZE must be +1
    void optimizeMap8(uint8_t *MAPA , uint8_t  *MAPAO, int lRef, int *LenOpMap){
    	unsigned int posRel=1, posVec=0;
    	uint64_t mask=0x01, aux=0x0;
    	/*Generando el mapa optimizado*/
    	for  (int y=1;y<lRef;y++){  //ojo diferencias reerse y forwad

    		if (aux!=0) aux=aux<<1;
    		if(MAPA[y]==1){ //Meta un 1
    			aux=mask|aux;

    		}

    		if ((posRel>=8)||(y>=lRef)) {
    			if (y>=lRef) aux=aux<<(8-posRel); //o -1?
    			MAPAO[posVec]=aux;
    			//if (y<=9) printf ("  iter  %i AUX %"PRIu64 "\n",y, aux);
    			//if (y<=9) printf ("  iter  %i AUX %lu\n",y, aux);//both are the same

    			if (y<lRef){
    				 posRel=1; 			posVec++;   			aux=0x0;

    			}
    		} else {
    			posRel++;
    		}

    	}//For Int y

    	*LenOpMap=posVec;
    } /*Fin del mapa optimizado*/



  char* itoa(int val, int base){

  	static char buf[32] = {0};

  	int i = 30;

  	for(; val && i ; --i, val /= base)

  		buf[i] = "0123456789abcdef"[val % base];

  	return &buf[i+1];

  }


  /*
//Decodificar: Convierte un Cromosoma x en su equivalente numero real
double Decodificar(Cromosoma x) {return ((double)x/(16777215.0));}
//fin Decodificar

//ConvertCrombin: Convierte Cromosoma entero a binario
//a: cromosoma a convertir
//b: cadena de numeros binarios (salida)
//c: cantidad de bits a usar
void ConvertCromBin (Cromosoma a, char b[], int c){
	int u;
	for (u=0; u<c-1; u++) b[u]='0';
   b[u]='\0';
   u=c-2;
	while (a!=0){
   	if((a%2)==1) b[u]='1';
      a /= 2;
      u--;
   }//fin mientras a != 0
}//fin Conversion Cromosoma -> Binario

revisar si el transcriptor de mapa orptimizado llega a la cantidad correcta de POSICIONES
Determinar Deleciones contiguas (hasta 4)
inserciones repetidas contiguas. Cuidado en relacion a la letra exacta del read.
susticituiones repetidas contiguas.  Cuidado en relacion a la letra exacta del read.
 *
 * */


  /*
   * g++ -std=c++11 -g3 -I ~/include -L ~/lib fm-index.cpp aligner.cpp extend.cpp -o fm-index -lsdsl -ldivsufsort -ldivsufsort64
   *
   * */
