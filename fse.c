#include <getopt.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>





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
	long dnBits;
};


/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
void compress_image(char* input_filename, FILE* output_file, long q) {

long *qmap_small = 0;
long *qmap = 0;


long** image=0;
long** quantised=0;

long nx, ny;
long i,j;
long* in_buff; 

  /* read image into buffer */
read_pgm_and_allocate_memory(input_filename, &nx, &ny, &image);

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
 printf("here \n");
//converting 2D to 1D
 long k=0;
 	for(i=1;i<=nx;i++)
 		for(j=1;j<=ny;j++)
 			{
 				in_buff[k]=image[i][j];
 				k++;
 			}
 	long srcSize=nx*ny;
//finding optimal value of tableLog
 	//minTableLog should be chosen from minBitsSrc and minBitsSymbols, as the source size > no.of symbols
 	// skip the computation and chose minBitsSrc as minBits, but 11>log(255)+2, and maxBitsSrc>11, we choose 11 as default tableLog
 	// by making some assumptions about srcSize and no. of symbols, the default value in the original code is 11, and here we skip
 	// the computations and just choose 11  (refer fse_compress -> optimalTableLog)


long *count,*norm_count,total=0,total_norm=0,max=0,tableSize=2048,maxValue=256,tableLog=log2(tableSize);
alloc_vector_int(&count,256);
alloc_vector_int(&norm_count,256);
for(i=0;i<maxValue;i++)
	count[i]=0;
for(i=0;i<srcSize;i++)
	count[in_buff[i]]++;
for(i=0;i<maxValue;i++)
{
	total+=count[i];
	printf("\n %ld -> %ld ", i, count[i]);
}
double t;
for(i=0;i<maxValue;i++)
{
	
	
	t=((double)count[i]/(double)total)*(double)tableSize;
	norm_count[i]=(fabs(t-floor(t)<0.5))?(long)floor(t):(long)ceil(t);
	total_norm+=norm_count[i];
	if(norm_count[i]>norm_count[max])
		max=i;

}
norm_count[max]+=(tableSize-total_norm); //correcting for rounding errors

total_norm=0;
for(i=0;i<maxValue;i++)
{
	printf("\n %ld -> %ld ", count[i], norm_count[i]);
	total_norm+=norm_count[i];
}
printf("\n M=%ld N=%ld", total_norm,total);

//************building the table****************

//symbol start positions 
long *cumul;
alloc_vector_int(&cumul,maxValue+2);
cumul[0];
for(i=1;i<maxValue;i++)
	cumul[i]=cumul[i-1]+norm_count[i-1];
cumul[maxValue]=tableSize+1;


printf("symbol start positions \n");
    for(i=0;i<=maxValue;i++)
    	printf("%lu ",cumul[i]);



 	
//#define FSE_TABLESTEP(tableSize) ((tableSize>>1) + (tableSize>>3) + 3) in original code
   // used to define the step size to distribute symbols given the table size. Have no idea why???
//substitute 11 and get step=1283
//now, we distribute the symbols
long *tableSymbol,step, pos=0;
step=(tableSize/2)+(tableSize/8)+3;
alloc_vector_int(&tableSymbol,tableSize);
for(i=0;i<256;i++) //iterating over symbols
{
	for(j=0;j<norm_count[i];j++)  //iterating for the number of normalized occurences
	{
		tableSymbol[pos]=i;
		pos=(pos+step)%tableSize;
	}
 }

 for(i=0;i<2048;i++)
    	printf("\n %ld -> %ld",i,tableSymbol[i] );

struct symbolTransform symbolTT[maxValue];
long maxBitsOut,minStatePlus;
total=0;
for(i=0;i<255;i++)
{
	switch(norm_count[i])
	{
		case 0:
		symbolTT[i].dnBits=(tableLog+1)*pow(2,16)-pow(2,tableLog);
		symbolTT[i].dstate=0;
		break;

		case 1:
		symbolTT[i].dnBits=(tableLog+1)*pow(2,16)-pow(2,tableLog);
		symbolTT[i].dstate=total-1;
		total++;
		break;

		default:
		maxBitsOut = (long)(tableLog - log(norm_count[i]-1));
        minStatePlus = (long)(norm_count[i]*pow(2, maxBitsOut));
        symbolTT[i].dnBits = (maxBitsOut *pow(2, 16)) - minStatePlus;
        symbolTT[i].dstate = total - norm_count[i];
        total +=  norm_count[i];


	}
}

printf("\n symbol transformation\n");
    for (i=0; i<maxValue; i++)
    	printf("%ld->%ld,%ld \n", i,symbolTT[i].dstate,symbolTT[i].dnBits);

 

  
}


/*--------------------------------------------------------------------------*/
void decompress_image(FILE* input_file, char* output_filename) {
  
}
  
/*--------------------------------------------------------------------------*/

int main(int argc, char** args) {
  int ch;	  /* all-purpose char */
  char used[256]; /* flag for all possible chars, used to check for multiple
                   * occurences of input parameters */

  char* input_filename=0;
  char* compressed_filename=0;
  char* output_filename=0;

  FILE* input_file=0;
  FILE* output_file=0;
  long q=256;



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
  
  output_file = fopen(compressed_filename, "wb");
  compress_image(input_filename,output_file,q);
  
  //fclose(input_file);
  fclose(output_file);

/*  input_file = fopen(compressed_filename, "rb");
  decompress_image(input_file,output_filename);
  fclose(input_file);
  */
  return 0;
}
