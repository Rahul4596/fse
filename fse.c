#include <getopt.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "lib/bfio.h"





/*----------------------------------------------------------------------------*/
float entropy
(
 long** u,                 /* image */
 long nx,                    /* x dimension */
 long ny                     /* y dimension */
 )
/* computes entropy of image */
{
  long i,j;

  float entropy  = 0.0f;
  long h[256];
  
  float p;
  for (i=0;i<=255;i++) {
    h[i]=0;
  }

  /* compute histogram */
    for (i=1; i<=nx; i++) {
      for (j=1; j<=ny; j++) {
        h[u[i][j]]++;
      }
    }
  
  
  for (i=0;i<=255;i++) {
    if (h[i]>0) {
      printf("%ld: %ld\n",i,h[i]);
      p = (float)h[i]/(float)(nx*ny);
      entropy -= p*logf(p)/logf(2.0);
    }
  }

  return entropy;
} /* entropy */




/*--------------------------------------------------------------------------*/

void alloc_vector_int

     (long **vector,   /* vector */
      long  n1)         /* size */

     /* allocates memory for a vector of size n1 */


{
*vector = (long *) malloc (n1 * sizeof(long));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough memory available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (long ***matrix,  /* matrix */
      long  n1,         /* size in direction 1 */
      long  n2)         /* size in direction 2 */
     /* allocates memory for matrix of size n1 * n2 */
{
long i;

*matrix = (long **) malloc (n1 * sizeof(long *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough memory available\n");
   exit(1);
   }
for (i=0; i<n1; i++)
    {
    (*matrix)[i] = (long *) malloc (n2 * sizeof(long));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough memory available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory

     (const char  *file_name,    /* name of pgm file */ 
      long        *nx,           /* image size in x direction, output */
      long        *ny,           /* image size in y direction, output */
      long        ***u)          /* image, output */   

/* 
  reads a greyscale image that has been encoded in pgm format P5;
  allocates memory for the image u; 
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
FILE   *inimage;    /* input file */
char   row[80];     /* for reading data */
long   i, j;        /* loop variables */

/* open file */
inimage = fopen (file_name, "rb");
if (NULL == inimage) 
   {
   printf ("could not open file '%s' for reading, aborting.\n", file_name);
   exit (1);
   }

/* read header */
fgets (row, 80, inimage);          /* skip format definition */
fgets (row, 80, inimage);        
while (row[0]=='#')                /* skip comments */
      fgets (row, 80, inimage);
sscanf (row, "%ld %ld", nx, ny);   /* read image size */
fgets (row, 80, inimage);          /* read maximum grey value */

/* allocate memory */
alloc_matrix (u, (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=1; j<=(*ny); j++) 
 for (i=1; i<=(*nx); i++) 
     (*u)[i][j] = (long) getc(inimage);

/* close file */
fclose(inimage);

return;

} /* read_pgm_and_allocate_memory */

/*--------------------------------------------------------------------------*/


void write_pgm

     (long  **u,          /* image, unchanged */ 
      long   nx,           /* image size in x direction */
      long   ny,           /* image size in y direction */
      char   *file_name,   /* name of pgm file */
      char   *comments)    /* comment string (set 0 for no comments) */

/* 
  writes a greyscale image into a pgm P5 file;
*/

{
FILE           *outimage;  /* output file */
long           i, j;       /* loop variables */
float          aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage) 
   {
   printf("Could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
   }

/* write header */
fprintf (outimage, "P5\n");                  /* format */
if (comments != 0)
   fprintf (outimage, comments);             /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     aux = u[i][j] + 0.499999;    /* for correct rounding */
     if (aux < 0)
        byte = (unsigned char)(0.0);
     else if (aux > 255)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(aux);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }

/* close file */
fclose (outimage);

return;

} /* write_pgm */


/*--------------------------------------------------------------------------*/
void prepare_quantisation_map(long q, /* number of quantised values */
                              long max_q, /* maximum possible grey value */
                              long* quantisation_map_small,
                              long* quantisation_map_back) {
  long i;
  long a = 0;
  long max_levels=65536; /* 16 bit */

  //printf("preparing quantisation map with q=%ld and max_q=%ld\n",q,max_q); 
        
  if (q > max_levels) {
    printf("warning, quantisation exceeds max levels, set to %ld\n",max_levels);
    q = max_levels;
  }
        
  if ((q == max_levels || q == 0) && (max_q == max_levels)) {
      for (i = 0; i < q; ++i) {
        quantisation_map_back[i] = i;
        quantisation_map_small[i] = i;
      }
  }
  
  for (i = 0; i <= q; ++i) {
    quantisation_map_back[i] =
      (long)((i+0.5)*max_q/q);
  }
  //printf("map back: done\n");
        
  a=0;
  for (i = 0; i <= max_q-1; ++i) {
    if (i-quantisation_map_back[a] > quantisation_map_back[a+1]-i
        && a < q-1) a++;
    quantisation_map_small[i] = a;
  }
  //printf("map small: done\n");
}


/*--------------------------------------------------------------------------*/
void quantise_image_scalar(long** original,
                           long nx, long ny,
                           long max_q,
                           long q,
                           long* quantisation_map_small,
                           long* quantisation_map_back,
                           long** quantised) {
  
  long x, y;
  long bx = 1, by = 1;
  
  prepare_quantisation_map(q,max_q,quantisation_map_small,quantisation_map_back);

    for (x = bx; x < nx + bx; x++)
      for (y = by; y < ny + by; y++) {
        if (original[x][y]>255) {
          original[x][y]=255;
        }
        if (original[x][y]<0) {
          original[x][y]=0;
        }
        quantised[x][y] = quantisation_map_small[(long)(original[x][y])];
        /*if (quantised[c][x][y]>=q) {
          printf("%f > %ld\n",original[c][x][y],quantised[c][x][y]);
          }*/
      }
}


/*--------------------------------------------------------------------------*/
void init_image_quantised(long** quantised,
                          float** mask,
                          long nx, long ny,
                          long* qmap,
                          long** image) {
  
  long i, j;

  if (mask==0) {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
          image[i][j]=(float)qmap[quantised[i][j]];
      }
    }
  } else {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
        if (mask[i][j]>0.5) {
          image[i][j]=(float)qmap[quantised[i][j]];
        } else {
          image[i][j]=0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------------------*/
struct symbolTransform
{
	unsigned dstate;
	unsigned dBits;
};


/*--------------------------------------------------------------------------*/
void write_long_bitwise(long c, /* (positive) number to write */
                        long n, /* number of bits */
                        BFILE* output_file) /* 0 no output,
                                               otherwise write to binary file */
{
  
  while(n>0)
  {
    /*set_bit(output_file,(c>>(n-1))&1);
    printf("b=%ld  " , (c>>(n-1))&1);
  	n--;*/
  	set_bit(output_file,c%2);
  	//printf("b=%ld ", c%2);
  	c/=2;
  	n--;
  }
  return;
}

/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
void compress_image(char* input_filename, BFILE* output_file, long q, unsigned* norm_count, unsigned* finalState,long* nx_r,long* ny_r) {

long *qmap_small = 0;
long *qmap = 0;


long** image=0;
long** quantised=0;

long nx, ny;
long i,j;
long* in_buff; 

  /* read image into buffer */
read_pgm_and_allocate_memory(input_filename, &nx, &ny, &image);
*nx_r=nx;*ny_r=ny;
  /* quantise if desired */
  if (q < 256) {
    alloc_matrix (&quantised, nx+2, ny+2);
    alloc_vector_int(&qmap,257);
    alloc_vector_int(&qmap_small,257);
    quantise_image_scalar(image,nx,ny,256,q,qmap_small,qmap,quantised);
    init_image_quantised(quantised,0,nx,ny,qmap,image);
  }
  alloc_vector_int(&in_buff,(nx*ny));
 //long minSize=5,maxSize=15;
 
//converting 2D to 1D
 long k=0;
 	for(i=1;i<=nx;i++)
 		for(j=1;j<=ny;j++)
 			{
 				in_buff[(i-1)+(j-1)*nx]=image[i][j];
 				//printf("(%ld,%ld)->%ld \n ", i,j, image[i][j]);
 				k++;
 			}
 	long srcSize=nx*ny;

 	//printf("here %ld \n", in_buff[k-1]);
//finding optimal value of tableLog
 	//minTableLog should be chosen from minBitsSrc and minBitsSymbols, as the source size > no.of symbols
 	// skip the computation and chose minBitsSrc as minBits, but 11>log(255)+2, and maxBitsSrc>11, we choose 11 as default tableLog
 	// by making some assumptions about srcSize and no. of symbols, the default value in the original code is 11, and here we skip
 	// the computations and just choose 11  (refer fse_compress -> optimalTableLog)


long *count,total=0,total_norm=0,max=0,tableSize=2048,maxValue=256,tableLog=log2(tableSize);
alloc_vector_int(&count,256);
unsigned temp_norm[maxValue];
//alloc_vector_int(&norm_count,256);
for(i=0;i<maxValue;i++)
	count[i]=0;
for(i=0;i<srcSize;i++)
	count[in_buff[i]]++;
for(i=0;i<maxValue;i++)
{
	total+=count[i];
	//printf("\n %ld -> %ld ", i, count[i]);
}
double t;
for(i=0;i<maxValue;i++)
{
	
	
	t=((double)count[i]/(double)total)*(double)tableSize;
	norm_count[i]=(long)floor(t);//(fabs(t-floor(t)<0.5))?(long)floor(t):(long)ceil(t);
  if(norm_count[i]==0)
    norm_count[i]=1;
	total_norm+=norm_count[i];
	if(norm_count[i]>norm_count[max])
		max=i;

}

printf("%ld %ld \n", norm_count[max],total_norm);
/*norm_count[max]+=(tableSize-total_norm); //correcting for rounding errors

total_norm=0;
for(i=0;i<maxValue;i++)
{
	printf("\n %ld -> %ld ", i, norm_count[i]);
  temp_norm[i]=norm_count[i];
	total_norm+=norm_count[i];
}*/


////////////////////////////////////
/*norm_count[20]++;
norm_count[36]--;
norm_count[44]--;
norm_count[52]--;
norm_count[92]--;
//norm_count[100]=181;
norm_count[124]--;
norm_count[140]--;
norm_count[172]--;
norm_count[180]--;
norm_count[204]--;
norm_count[220]--;
norm_count[252]--;
total_norm-=10;*/

///////////////////////////////////

int flag=0;
if(tableSize-total_norm>0)
  norm_count[max]+=(tableSize-total_norm); //correcting for rounding errors
else
{
  while(1)
  {
    for(i=0;i<256;i++)
    {
      if(norm_count[i]>1)
        {
          norm_count[i]--;
          total_norm--;
          printf("%ld \n",total_norm-tableSize );
          if(total_norm-tableSize==0)
            {flag=1;break;}
        }
    }
    if(flag)
      break;
  }
}
total_norm=0;
for(i=0;i<maxValue;i++)
{
  if(count[i])
  printf("\n %ld -> %ld ", count[i], norm_count[i]);
  temp_norm[i]=norm_count[i];
  total_norm+=norm_count[i];
}

///////////////////////////////////////////////////////////

printf("\n M=%ld N=%ld", total_norm,total);

//************building the table****************

//symbol start positions 
long *cumul;
alloc_vector_int(&cumul,maxValue+2);
cumul[0];
for(i=1;i<maxValue;i++)
	cumul[i]=cumul[i-1]+norm_count[i-1];
cumul[maxValue]=tableSize+1;


//printf("symbol start positions \n");
//    for(i=0;i<=maxValue;i++)
//    	printf("%lu ",cumul[i]);



 	
   // used to define the step size to distribute symbols given the table size. 
    //Have some idea now. As symbols should be spread as far as possible from each other
    //ie (0,1/2,1/4,3/4,....), any sequence with step size n/2+n/4+n/8+...+odd_number
    // will approimate this distribution. Also, as the step size and n value are coprime
    // all states will be visited exactly once.

 //#define FSE_TABLESTEP(tableSize) ((tableSize>>1) + (tableSize>>3) + 3) in original code
//substitute 11 and get step=1283
//now, we distribute the symbols
long step, pos=0;
step=(tableSize/2)+(tableSize/8)+3;
long *tableSymbol;
alloc_vector_int(&tableSymbol,tableSize);
for(i=0;i<256;i++) //iterating over symbols
{
	for(j=0;j<norm_count[i];j++)  //iterating for the number of normalized occurences
	{
		tableSymbol[pos]=i;
		pos=(pos+step)%tableSize;
	}
 }


 
 /*for(i=0;i<tableSize;i++)
  {
    unsigned sym=tableSymbol[i];
    
    unsigned nextState= temp_norm[sym]++;
    unsigned nBits = (unsigned)(tableLog - log2(nextState)+1);
    if(nBits>tableLog)
      nBits=tableLog;
    //if(i==357)
    //  decode_table[i].nBits--;
    unsigned newState=(unsigned)((nextState<< nBits) - tableSize);
    if(newState==tableSize)
    {
      /*if(temp_norm[sym]==1)
        {norm_count[sym]++;total_norm++;}
      else
      if(norm_count[sym]!=1)
        {norm_count[sym]--;total_norm--;}

    }
    printf("state->%ld symbol->%ld nBits->%ld nextState->%ld newState->%ld \n", i, sym, nBits, nextState, newState);
    
  }

  norm_count[max]+=(tableSize-total_norm);
  total_norm+=(tableSize-total_norm);

  for(i=0;i<tableSize;i++)
  {
    unsigned sym=tableSymbol[i];
    
    unsigned nextState= temp_norm[sym]++;
    unsigned nBits = (unsigned)(tableLog - log2(nextState)+1);
    if(nBits>tableLog)
      nBits=tableLog;
    //if(i==357)
    //  decode_table[i].nBits--;
    unsigned newState=(unsigned)((nextState<< nBits) - tableSize);
printf("state->%ld symbol->%ld nBits->%ld nextState->%ld newState->%ld \n", i, sym, nBits, nextState, newState);
    
  }
  

 /*for(i=0;i<256;i++)
    	{
        total_norm+=norm_count[i];
      }*/
   // printf("\n M=%ld N=%ld", total_norm,total);


// we build the table. The table is an arrya of struct.Each symbol has 3 fields
 //minBits, minState-> if the state value is higher than this, use minBits+1
 // dstate gives the change in state after encoding

 /* Build table */
unsigned tableU16[4200];
 
    {   long u; for (u=0; u<tableSize; u++) {
        long s = tableSymbol[u];   /* note : static analyzer may not understand tableSymbol is properly initialized */
        tableU16[cumul[s]++] = (unsigned) (tableSize+u);   /* TableU16 : sorted by symbol order; gives next state value */
    }   }

struct symbolTransform symbolTT[maxValue];
long maxBitsOut,minStatePlus;
total=0;
for(i=0;i<255;i++)
{
	switch(norm_count[i])
	{
		case 0:
		//symbolTT[i].dnBits=(tableLog+1)*pow(2,16)-pow(2,tableLog);
		symbolTT[i].dstate=0;
		break;
		//case -1:
		case 1:
		symbolTT[i].dBits=(tableLog<<16) - (1<< (tableLog-1));
		//symbolTT[i].minStat=pow(2,tableLog);
		symbolTT[i].dstate=total-1;
		//printf("\n symbol=%ld  dstate=%ld dBits=%ld",i,symbolTT[i].dstate, symbolTT[i].dBits);
	//	if(i==20){
     

      //printf("symbol= %ld",i);
      //int s;
      //for(s=0;s<2048;s++)
      //printf("%ld->%ld ", s,(symbolTT[i].minBits+s))>>16;
	//	}
		total++;
		break;

		default: ;
		unsigned maxBits = (long)(tableLog - (long)log2(norm_count[i]-1));
    //if(maxBits>tableLog) maxBits=tableLog;
        unsigned minStatePlus = (long)(norm_count[i]<< (maxBits));
        //symbolTT[i].dnBits = (maxBitsOut *pow(2, 16)) - minStatePlus;
        symbolTT[i].dstate = total - norm_count[i];
        symbolTT[i].dBits= (maxBits << 16) - minStatePlus;
        total +=  norm_count[i];
        printf("\n symbol=%ld maxBitsOut=%ld minState=%ld dstate=%ld dBits=%ld",i,maxBits, minStatePlus,symbolTT[i].dstate, symbolTT[i].dBits);
       /* if(i==92)
    {
    */  
   /* }*/

	}
	//if(norm_count[i])
	//printf("\n symbol=%ld minBitsOut=%ld minStatePlus=%ld dstate=%ld",i,symbolTT[i].minBits,symbolTT[i].minStatePlus, symbolTT[i].dstate);
}

/*int s;
      for(s=1;s<2048;s++)
      printf("%ld->%ld \n",s,(long)log2(s));*/

//printf("k=%ld \n",k);
struct symbolTransform symbolT;
unsigned stateValue=tableSize;
unsigned nBitsOut;
unsigned out_buff[nx*ny];


for(i=k-1;i>=0;i--)
{
	symbolT=symbolTT[in_buff[i]];
	//nBitsOut=stateValue<symbolT.minStatePlus?symbolT.minBits:symbolT.minBits+1;
	if(i<nx*ny-1)
		out_buff[i]=stateValue;
	nBitsOut = (unsigned)(stateValue+symbolT.dBits) >> 16;
	

	//printf("i=%ld symbol=%ld %ld,%ld \n",i,in_buff[i],stateValue,nBitsOut );
	//encode state value push bits(state value)
	//write_long_bitwise(stateValue,nBitsOut,output_file);
	
	stateValue= tableU16[(long)((stateValue>>(nBitsOut))+symbolT.dstate)];
	if(i>=nx*ny-500)
	{printf("k=%ld symbol=%ld %ld,%ld \n",i,in_buff[i],stateValue,nBitsOut );}

}
out_buff[nx*ny-1]=stateValue;
*finalState=stateValue;

for(i=nx*ny-1;i>=0;i--)
	{
		
		symbolT=symbolTT[in_buff[nx*ny-1-i]];
		nBitsOut = (unsigned)(out_buff[nx*ny-1-i]+symbolT.dBits) >> 16;
		write_long_bitwise(out_buff[nx*ny-1-i],nBitsOut,output_file);
		
		//if(i>nx*ny-1-300)
			//printf("i=%ld symbol=%ld %ld,%ld \n",nx*ny-1-i, in_buff[nx*ny-1-i], out_buff[nx*ny-1-i],nBitsOut);
	}

/*printf("\n symbol transformation\n");
    for (i=0; i<maxValue; i++)
    	printf("%ld->%ld,%ld \n", i,symbolTT[i].dstate,symbolTT[i].dnBits);
*/
// start encoding from the last symbol and go to the first symbol , but first we have to
//initialise state


  
}



/*--------------------------------------------------------------------------*/
void decompress_image(FILE* input_file, char* output_filename, unsigned* norm_count, unsigned initialState, long nx, long ny)
{
	

	// first step in decoding is building the decoding table. For now, we assume we have the norm_count values
	// but later we have to write it in the file and read from it.

	//first spread the symbols

	long tableLog=11;
	long tableSize=pow(2,tableLog);
	long out_buff[nx*ny];
	struct decode_symbol
	{
		unsigned newState;
		unsigned symbol;
		unsigned nBits;
	} decode_table[tableSize];

	long pos=0;
	long step=(tableSize/2)+(tableSize/8)+3;
	long *tableSymbol;
	alloc_vector_int(&tableSymbol,tableSize);
	unsigned i,j;
	for(i=0;i<256;i++) //iterating over symbols
	{
		for(j=0;j<norm_count[i];j++)  //iterating for the number of normalized occurences
		{
			decode_table[pos].symbol=i;
			pos=(pos+step)%tableSize;
		}
	}

	//build decoding table

	for(i=0;i<tableSize;i++)
	{
		unsigned sym=decode_table[i].symbol;
		unsigned nextState= norm_count[sym]++;
		decode_table[i].nBits = (unsigned)(tableLog - (long)log2(nextState));
		if(decode_table[i].nBits>tableLog)
      decode_table[i].nBits=tableLog;
    
		decode_table[i].newState=(unsigned)((nextState<< decode_table[i].nBits) - tableSize);
    if(decode_table[i].newState==tableSize)
      {decode_table[i].newState=(unsigned)(((nextState-1)<< decode_table[i].nBits) - tableSize); }
		//printf("state->%ld symbol->%ld nBits->%ld nextState->%ld newState->%ld \n", i, decode_table[i].symbol, decode_table[i].nBits, nextState, decode_table[i].newState);
		
	}

	unsigned stateValue=initialState-tableSize;
	//printf("%ld \n", stateValue);
	struct decode_symbol stateInfo;
	long nextBits;
	long k;
	for(k=0;k<nx*ny;k++)
	{
			stateInfo=decode_table[stateValue];
		out_buff[k]=stateInfo.symbol;

    if(k>=nx*ny-500)
    printf("%ld->%ld %ld %ld",k,out_buff[k], stateValue,stateInfo.nBits );
		nextBits=get_bits(input_file,stateInfo.nBits);
		
		stateValue=(stateInfo.newState + nextBits) ;
		 

	}

	long** image;
	alloc_matrix (&image, nx+2, ny+2);
  
  	for (i=1;i<=nx;i++) {
    	for (j=1;j<=ny;j++) {
      image[i][j]=(long)out_buff[(i-1)+(j-1)*nx];
      
    }
  }

  //entropy(image,nx,ny);

  
  write_pgm(image,nx,ny,output_filename,0);


	// */



}
  

  
/*--------------------------------------------------------------------------*/

int main(int argc, char** args) {
  int ch;	  /* all-purpose char */
  char used[256]; /* flag for all possible chars, used to check for multiple
                   * occurences of input parameters */

  char* input_filename=0;
  char* compressed_filename=0;
  char* output_filename=0;

  //FILE* input_file=0;
  //FILE* output_file=0;
  long q=256;

  BFILE* output_file=0;

  BFILE* input_file=0;
long nx,ny;


  /* process input parameters */
  for (ch = 0;ch <= 255; ch++)
    used[ch] = 0;
  while ((ch = getopt(argc,args,"i:o:c:q:")) != -1) {
    used[ch]++;

    if (used[ch] > 1) {
      printf("Double parameter: %c\n",ch);
      printf("Please check your input again.\n");
      
    }

    switch(ch) {
      case 'i': input_filename = optarg;break;
      case 'c': compressed_filename = optarg;break;
      case 'o': output_filename = optarg;break;
      case 'q': q = atoi(optarg);break;
      default: printf("Unknown argument.\n");
    }
  }

  if (input_filename==0) {
    printf("Image file must be specified\n");
    return 0;
  }

  /*input_file = fopen(input_filename, "rb");*/
   
  unsigned finalState;
  unsigned norm_count[256];
  long tableSymbol[2048];
  output_file = bfopen(compressed_filename, "w");
  compress_image(input_filename,output_file,q,norm_count,&finalState,&nx,&ny);
  //write_long_bitwise(2693,4,output_file);
  //fclose(input_file);
  bfclose(output_file);
  
  input_file = bfopen(compressed_filename, "r");
  decompress_image(input_file,output_filename, norm_count,finalState,nx,ny);
  bfclose(input_file);


  
  return 0;
}
