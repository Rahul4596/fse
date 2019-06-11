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

unsigned mylog2(long x)
{
  unsigned n=0;
  while(x)
  {
    x=x>>1;
    n++;
  }
  return n-1;
}

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

typedef struct _BITBUFFER {
  unsigned char *b; /* array for bytewise content */
  long max_size;    /* maximum size of the buffer (determined by allocation) */
  long byte_pos;    /* current byte position */
  long bit_pos;     /* current bit position (inside current byte) */
} BITBUFFER;

void alloc_bitbuffer

  (BITBUFFER *buffer, /* bitbuffer */
   long n)            /* size */

/* allocates memory for a bitbuffer of size n */
{
long i; /* loop variable */

/* allocate array for individual bytes */
buffer->b = (unsigned char *) malloc(n * sizeof(unsigned char));
if (buffer->b == NULL)
{
   printf("alloc_bitbuffer: not enough memory available\n");
   exit(1);
}

/* set max size, positions, and initialise all bytes with 0 */
buffer->max_size = n;
for (i = 0; i < buffer->max_size; i++)
{
   buffer->b[i] = 0;
}
buffer->byte_pos = 0;
buffer->bit_pos = 0;
return;
} /* alloc_bitbuffer */

/*--------------------------------------------------------------------------*/

long bitbuffer_addbit

  (BITBUFFER *buffer,  /* bitbuffer (output) */
   unsigned char bit)  /* bit to be written */

/* Add a single bit to the buffer at the current bit position.
 * Returns 1 on success, 0 if the buffer is full
 */
{
if ((buffer->byte_pos < buffer->max_size) && (buffer->bit_pos < 8))
{
   if (bit)
      buffer->b[buffer->byte_pos] |= (1 << (buffer->bit_pos));
   else
      buffer->b[buffer->byte_pos] &= ~(1 << (buffer->bit_pos));
   (buffer->bit_pos)++;
   if (buffer->bit_pos > 7)
   {
      buffer->bit_pos = 0;
      buffer->byte_pos++;
   }
   return 1;
}
return 0; /* adding not successful, memory full */
} /* bitbuffer_addbit */

void write_long_bitwise_2(long c, /* (positive) number to write */
                        long n, /* number of bits */
                        BITBUFFER *buffer) /* 0 no output,
                                               otherwise write to binary file */
{
  
  
 //printf("writing->%ld with %ld bits \n",c,n );
  unsigned char bit;
  while(n--)
  {
   
    bit=(char)((c>>n)&0x1);
  if ((buffer->byte_pos < buffer->max_size) && (buffer->bit_pos < 8))
  {
   if (bit)
      buffer->b[buffer->byte_pos] |= (1 << (buffer->bit_pos));
   else
      buffer->b[buffer->byte_pos] &= ~(1 << (buffer->bit_pos));


   (buffer->bit_pos) = (buffer->bit_pos + 1)& 0x7;

  buffer->byte_pos+=(buffer->bit_pos == 0x0);
   
  }

    //printf("%ld ", (c>>n)&0x1);
    //c/=2;
    //n--;
  }
  
}

/*--------------------------------------------------------------------------*/

long bitbuffer_getbit

  (BITBUFFER *buffer)  /* bitbuffer (output) */

/* Get a single bit from the buffer at the current bit position and
 * move the position one bit further.
 * Returns the bit on success, -1 on failure.
 */
{
unsigned char bit;

if ((buffer->byte_pos < buffer->max_size) && (buffer->bit_pos < 8))
{
   bit = (buffer->b[buffer->byte_pos] >> buffer->bit_pos) & 1;
   (buffer->bit_pos)++;
   if (buffer->bit_pos > 7)
   {
      buffer->bit_pos = 0;
      buffer->byte_pos++;
   }
   return bit;
}
else
{
   return -1; /* end of buffer already reached, no more bits to read */
}
} /* bitbuffer_getbit */


long bitbuffer_getbits

  (BITBUFFER *buffer, long n)  /* bitbuffer (output) */

/* Get a single bit from the buffer at the current bit position and
 * move the position one bit further.
 * Returns the bit on success, -1 on failure.
 */
{

long num=0;
/*unsigned int temp=( (unsigned int)(buffer->b[buffer->byte_pos])<<16) +( (unsigned int)(buffer->b[buffer->byte_pos+1]) << 8) + (unsigned int)(buffer->b[buffer->byte_pos+2]);// + ((unsigned long)(buffer->b[((buffer->byte_pos) +1)])<<8) ;//+ (unsigned long)(buffer->b[(buffer->byte_pos + 2)]);


printf("%ld %ld %ld %lu %ld\n ",buffer->b[buffer->byte_pos],buffer->b[(buffer->byte_pos +1)], buffer->b[(buffer->byte_pos + 2)], temp,n);
temp = ((temp)<<(8 + buffer->bit_pos))>>(32-n);

long byt2= buffer->byte_pos +((buffer->bit_pos + n -1) >> 3);
long bit2= (buffer->bit_pos + n) & 0x7;*/

while(n--)
{
  
  num=(num<<1)+((buffer->b[buffer->byte_pos]>>buffer->bit_pos) &1);
  (buffer->bit_pos) = (buffer->bit_pos + 1)& 0x7;

  buffer->byte_pos+=(buffer->bit_pos == 0x0);

}


//printf("%lu,%lu,%lu,%lu,%lu,%lu \n",num,buffer->byte_pos,buffer->bit_pos,temp,byt2,bit2);


return num;
} /* bitbuffer_getbit */

/*--------------------------------------------------------------------------*/

void bitbuffer_writefile

  (BITBUFFER *buffer, /* buffer to be written */
   FILE *bitfile)     /* output file (should be open for writing) */

/* Write full bitbuffer to file in a single disk operation. */
{
fwrite(buffer->b, sizeof(unsigned char), buffer->byte_pos + 1, bitfile);
} /* bitbuffer_writefile */

/*--------------------------------------------------------------------------*/

void bitbuffer_loadfile

  (BITBUFFER *buffer, /* output buffer*/
   FILE *bitfile)     /* input file (should be open for reading "rb") */

/* Load file content from current position until end into bitbuffer. */
{
long start_pos; /* initial position of file pointer */
long chunksize; /* size of file chunk to read */


start_pos = ftell(bitfile);
fseek(bitfile, 0, SEEK_END);
chunksize = ftell(bitfile)-start_pos;
fseek(bitfile, start_pos, SEEK_SET);

/*printf("allocating buffer of size %ld\n",chunksize);*/
alloc_bitbuffer(buffer,chunksize);

/*printf("start reading at position %ld, write to %ld",
   start_pos,buffer->byte_pos); */

fread(buffer->b, 1, chunksize, bitfile);
/*printf("%ld bytes loaded\n",chunksize);*/

} /* bitbuffer_loadfile */





/* apply WNC algorithm for adaptive arithmetic integer encoding */
void encode_adaptive_wnc(
  long* sourceword,   /* array containing n numbers from {0,...,s}
                         where s is the end of file symbol */
  long n,             /* length of sourceword */
  long s,             /* size of source alphabet */
  float r,           /* rescaling parameter */
  long  M,            /* WNC discretisation parameter M */
  FILE* debug_file,   /* 0 - no output, otherwise debug output to file */
  FILE* compressed)  {/* binary file for compressed bitstring */

  long i,j;      /* loop variables */
  long L,H;      /* low and high interval endpoints of current interval */
  long oldL;     /* temporary variable to preserve L for computing new interval*/
  long C;        /* sum of all counters */
  long k;        /* underflow counter */
  long* counter; /* array of counters for adaptive probabilities */
  long M12=M/2, M14=M/4, M34=3*M/4;  /* time savers */
  long symbol;   /* index of current symbol in counter array */
  long csum;     /* sum of counters 0,...,symbol-1 */
  long rescaling_counter = 0;
  long bits_written = 0;


  /* allocate memory */
  alloc_vector_int(&counter,s);
  BITBUFFER buffer;
  alloc_bitbuffer(&buffer,16*n);
 // printf("%ld \n",8*n );


  /* initialise counters and C */
  for (i=0;i<s;i++) {
    counter[i]=1;
  }
  C = s;

  if (debug_file != 0) {
    fprintf(debug_file,"n: %ld, s: %ld, r: %f, M: %ld\n",n,s,r,M);
  }

  if ((double)C>(double)M/4.0+2.0) {
    if (debug_file != 0) {
      fprintf(debug_file,
              "M=%ld is too small (C=%ld), setting M to %ld\n",M,C,8*C);
    }
    M=8*C;
    M12=M/2; M14=M/4; M34=3*M/4;
  }

  /* initialise interval endpoints and underflow counter */
  L=0;
  H=M;
  k=0;

  /* encode sourceword */
  for (i=0;i<n;i++) {
    if (debug_file != 0) {
      fprintf(debug_file,"sourceword[%ld]=%ld\n",i,sourceword[i]);
    }
    /* underflow expansions/rescaling */
    while (1) {
      /* check for underflow expansion x -> 2*x - M/2 */
      if ((L >= M14) && (L<M12) && (H>M12) && (H<=M34)) {
        L=2*L-M12; H=2*H-M12;
        k++;
        if (debug_file != 0) {
          fprintf(debug_file,"underflow: x -> 2*x - %ld\n",M12);
          rescaling_counter++;
        }
        continue;
      }

      /* check for rescaling x -> 2*x */
      if (H<=M12) {
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x:");
          rescaling_counter++;
        }
        L=2*L; H=2*H;
        /* write 01^k to bitstream */
        bitbuffer_addbit(&buffer,(char)0);
        if (debug_file != 0) {
          fprintf(debug_file,"written bits: 0");
        }
        for (j=0;j<k;j++) {
          bitbuffer_addbit(&buffer,(char)1);
          if (debug_file != 0) {
            fprintf(debug_file,"1");
          }
        }
        if (debug_file != 0) {
          fprintf(debug_file,"\n");
        }
        k=0;
        continue;
      }

      /* check for rescaling x -> 2*x - M/2 */
      if (L>=M12) {
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x - %ld:",M);
          rescaling_counter++;
        }
        L=2*L-M; H=2*H-M;
        /* write 10^k to bitstream */
        bitbuffer_addbit(&buffer,(char)1);
        if (debug_file != 0) {
          fprintf(debug_file,"written bits: 1");
        }
        for (j=0;j<k;j++) {
          bitbuffer_addbit(&buffer,(char)0);
          if (debug_file != 0) {
            fprintf(debug_file,"0");
          }
        }
        if (debug_file != 0) {
          fprintf(debug_file,"\n");
        }
        k=0;
        continue;
      }

      if (debug_file != 0) {
        fprintf(debug_file,"k: %ld, [%ld, %ld)\n",k,L,H);
      }
      break;
    }

    /* readjustment */
    while ((double)C>(double)M/4.0+2.0) {
      if (debug_file != 0) {
        fprintf(debug_file,"C: %ld, M/4+2.0=%f\n",C,M/4.0+2.0);
      }
      C=0;
      for (j=0;j<s;j++) {
        counter[j]=(long)round((double)counter[j]*r);
        if (counter[j]==0) counter[j]=1; /* no zero-counters allowed */
        C+=counter[j];
      }
    }

    /* encode symbol */
    symbol=sourceword[i];
    csum=0;
    for (j=0;j<symbol;j++) csum+=counter[j];

    oldL=L;
    L=L+(long)floor((double)(csum*(H-L))/(double)C);

H=oldL+(long)floor((double)((csum+counter[symbol])*(H-oldL))/(double)C);
    if (debug_file != 0) {
      fprintf(debug_file,"new [L,H) = [%ld,%ld)\n",L,H);
    }
    counter[symbol]++; C++;
  }

      /* underflow expansions/rescaling */
    while (1) {
      /* check for underflow expansion x -> 2*x - M/2 */
      if ((L >= M14) && (L<M12) && (H>M12) && (H<=M34)) {
        L=2*L-M12; H=2*H-M12;
        k++;
        if (debug_file != 0) {
          fprintf(debug_file,"underflow: x -> 2*x - %ld\n",M12);
        }
        continue;
      }

      /* check for rescaling x -> 2*x */
      if (H<=M12) {
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x:");
        }
        L=2*L; H=2*H;
        /* write 01^k to bitstream */
        bitbuffer_addbit(&buffer,(char)0);
        bits_written++;
        if (debug_file != 0) {
          fprintf(debug_file,"written bits: 0");
        }
        for (j=0;j<k;j++) {
          bitbuffer_addbit(&buffer,(char)1);
          bits_written++;
          if (debug_file != 0) {
            fprintf(debug_file,"1");
          }
        }
        if (debug_file != 0) {
          fprintf(debug_file,"\n");
        }
        k=0;
        continue;
      }

      /* check for rescaling x -> 2*x - M/2 */
      if (L>=M12) {
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x - %ld:",M);
        }
        L=2*L-M; H=2*H-M;
        /* write 10^k to bitstream */
        bitbuffer_addbit(&buffer,(char)1);
        bits_written++;
        if (debug_file != 0) {
          fprintf(debug_file,"written bits: 1");
        }
        for (j=0;j<k;j++) {
          bitbuffer_addbit(&buffer,(char)0);
          bits_written++;
          if (debug_file != 0) {
            fprintf(debug_file,"0");
          }
        }
        if (debug_file != 0) {
          fprintf(debug_file,"\n");
        }
        k=0;
        continue;
      }

      if (debug_file != 0) {
        fprintf(debug_file,"k: %ld, [%ld, %ld)\n",k,L,H);
      }
      break;
    }


  /* last step */
  if (debug_file != 0) {
    fprintf(debug_file,"last interval - written bits:");
  }
  if (L<M14) {
      bitbuffer_addbit(&buffer,(char)0);
      bits_written++;
      if (debug_file != 0) {
        fprintf(debug_file,"0");
      }
    for (j=0;j<k+1;j++) {
      bitbuffer_addbit(&buffer,(char)1);
      bits_written++;
      if (debug_file != 0) {
        fprintf(debug_file,"1");
      }
    }
  } else {
    bitbuffer_addbit(&buffer,(char)1);
    bits_written++;
    if (debug_file != 0) {
      fprintf(debug_file,"1");
    }
    for (j=0;j<k+1;j++) {
      bitbuffer_addbit(&buffer,(char)0);
      bits_written++;
      if (debug_file != 0) {
        fprintf(debug_file,"0");
      }
    }
  }

  /* write additional bits that have to be read by decoder */
  if (debug_file != 0) {
    fprintf(debug_file,"\n padding bits:");
  }
  for (i=0;i<log2(M)-bits_written;i++) {
    bitbuffer_addbit(&buffer,(char)0);
    if (debug_file != 0) {
      fprintf(debug_file,"0");
    }
  }

  if (debug_file != 0) {
    fprintf(debug_file,"\n rescalings/underflow expanions: %ld, additional bits_written: %ld\n",
            rescaling_counter,bits_written);
  }

  /* allocate memory */
  bitbuffer_writefile(&buffer,compressed);

  //disalloc_long_vector(counter,s);
}



long log2long(long x) {
  int logx = 0;
  while (x >>= 1) ++logx;
  return logx;
}
/*--------------------------------------------------------------------------*/

/* apply WNC algorithm for adaptive arithmetic integer encoding */
void decode_adaptive_wnc(
    FILE* compressed,  /* binary file with compressed bitstring */
    long n,             /* length of sourceword */
    long s,             /* size of source alphabet */
    double r,           /* rescaling parameter */
    long  M,            /* WNC discretisation parameter M */
    FILE* debug_file,   /* 0 - no output, 1 - debug output to file */
    long* sourceword)  {/* array containing n numbers from {0,...,s}
                           where s is the end of file symbol */
  long i,j;      /* loop variables */
  long L,H;      /* low and high interval endpoints of current interval */
  long oldL;     /* temporary variable to preserve L for computing new interval*/
  long C;        /* sum of all counters */
  long* counter; /* array of counters for adaptive probabilities */
  long M12=M/2, M14=M/4, M34=3*M/4;  /* time savers */
  long symbol;   /* index of current symbol in counter array */
  long csum;     /* sum of counters 0,...,symbol-1 */
  long w;        /* variable for finding correct decoding inverval */
  long v;        /* partial dyadic fraction */
  long b;        /* auxiliary variable for reading individual bits */
  long N;        /* auxiliary variable for determining number of initial bits
                    for v */
 
  printf("n=%ld s=%ld M=%ld \n",n,s,M );
  /* allocate memory */
  alloc_vector_int(&counter,s);



  BITBUFFER buffer;
  alloc_bitbuffer(&buffer,16*n);
  bitbuffer_loadfile(&buffer,compressed);
  
  /* initialise counters and C */
  for (i=0;i<s;i++) {
    counter[i]=1;
  }
  C = s;

  if ((double)C>(double)M/4.0+2.0) {
    if (debug_file != 0) {
      fprintf(debug_file,
              "M=%ld is too small (C=%ld), setting M to %ld\n",M,C,8*C);
    }
    M=8*C;
    M12=M/2; M14=M/4; M34=3*M/4;
  }

  /* initialise interval endpoints*/
  L=0;
  H=M;

  /* read first bits of codeword to obtain initival v */
  N=log2long(M); /* assumes that M is a power of 2! */
  j=(long)pow(2,N-1);
  v=0;
  for (i=0;i<N;i++) {
   b = bitbuffer_getbit(&buffer);
   v += b*j;
   if (debug_file != 0) {
     fprintf(debug_file,
             "v: %ld, j: %ld, i: %ld, b: %ld\n",v,j,i,b);
   }
   j/=2;
  }
  if (debug_file != 0) {
    fprintf(debug_file,"initial v: %ld (%ld first bits from coded file, %f)\n",
           v,N,log((double)M)/log(2.0));
  }

  


  /* decode sourceword */
  for (i=0;i<n;i++) {
  
    /* underflow expansions/rescaling */
    while (1) {
      /* check for underflow expansion x -> 2*x - M/2 */
      if ((L >= M14) && (L<M12) && (H>M12) && (H<=M34)) {
        L=2*L-M12; H=2*H-M12; v=2*v-M12;
        /* shift in next bit */
        b=bitbuffer_getbit(&buffer);
        if (b!=EOF) {
          v+=b;
        }
        if (debug_file != 0) {
          fprintf(debug_file,
                 "underflow: x -> 2*x - %ld, [%ld,%ld), b %ld v %ld\n",
                 M12,L,H,b,v);
        }
        continue;
      }
    
      /* check for rescaling x -> 2*x */
      if (H<=M12) {
        L=2*L; H=2*H; v=2*v;
        /* shift in next bit */
        b=bitbuffer_getbit(&buffer);
        if (b!=EOF) {
          v+=b;
        }
        if (debug_file != 0) {
          fprintf(debug_file,"rescaling: x-> 2*x, [%ld, %ld), b %ld, v %ld\n",
                  L,H,b,v);
        }
        continue;
      }
    
      /* check for rescaling x -> 2*x - M/2 */
      if (L>=M12) {
        L=2*L-M; H=2*H-M; v=2*v-M;
        /* shift in next bit */
        b=bitbuffer_getbit(&buffer);
        if (b!=EOF) {
          v+=b;
        }
        if (debug_file != 0) {
          fprintf(debug_file,
                  "rescaling: x-> 2*x - %ld, [%ld,%ld), b %ld, v %ld\n",
                 M,L,H,b,v);
        }
        continue;
      }

      if (debug_file != 0) {
        fprintf(debug_file,"v: %ld, [%ld, %ld)\n",v,L,H);
      }
      break;
    }
    
    /* readjustment */
    while ((double)C>(double)M/4.0+2.0) {
      if (debug_file != 0) {
        fprintf(debug_file,"readjust C: %ld, M/4+2.0=%f\n",C,M/4.0+2.0);
      }
      C=0;
      for (j=0;j<s;j++) {
        counter[j]=(long)round((double)counter[j]*r);
        if (counter[j]==0) counter[j]=1; /* no zero-counters allowed */
        C+=counter[j];
      }
    }

    /* decode symbol */
    w=((v-L+1)*C-1)/(H-L);

    /* find correct interval */
    symbol=0;
    csum=0;
    while ((symbol < s-1) && ((csum > w) || (csum+counter[symbol])<=w)) {
      csum+=counter[symbol];
      symbol++;
    }
    sourceword[i]=symbol;
    oldL=L;
    L=L+(long)floor((double)(csum*(H-L))/(double)C);
    H=oldL+(long)floor((double)((csum+counter[symbol])*(H-oldL))/(double)C);
    counter[symbol]++; C++;
    if (debug_file != 0) {
      fprintf(debug_file,"[c_i,c_i-1) = [%ld %ld) ",csum,csum+counter[symbol]);
      fprintf(debug_file,"w: %ld symbol[%ld]: %ld, new [L,H)=[%ld,%ld)\n",
              w,i,symbol,L,H);
    }
  }

}

void FSE_compress(long* in_buff, long srcSize, BITBUFFER *buffer)
{
  long *count,total=0,total_norm=0,max=0,tableSize=2048,maxValue=256,tableLog=mylog2(tableSize);

long *norm_count,i,j;


alloc_vector_int(&norm_count,256);



alloc_vector_int(&count,256); // must change if to use for runlengths
long num_sym=0;

//unsigned temp_norm[maxValue];

for(i=0;i<maxValue;i++)
  count[i]=0;
for(i=0;i<srcSize;i++)
  count[in_buff[i]]++;
for(i=0;i<maxValue;i++)
{
  total+=count[i];
  if(count[i])
    num_sym++;
  //printf("\n %ld -> %ld ", i, count[i]);
}
double t[maxValue];

for(i=0;i<maxValue;i++)
{
  
  
  if(count[i]==0)
  {
    norm_count[i]=0;
    continue;
  }
  t[i]=((double)count[i]/(double)total)*(double)tableSize;
  norm_count[i]=(long)floor(t[i]);//(fabs(t-floor(t)<0.5))?(long)floor(t):(long)ceil(t);
  if(norm_count[i]==0)
    norm_count[i]=1;
  total_norm+=norm_count[i];
  if(norm_count[i]>norm_count[max])
    max=i;


}





//printf(" total - %ld \n",total_norm );

int flag=0;
if(tableSize-total_norm>0)
  {
    flag=0;
  while(tableSize-total_norm>0)
  {
    for(i=1;i<maxValue;i++)
    {
      //printf("%ld ",i );
      if(t[i]-floor(t[i]) > t[flag] - floor(t[flag]) || (t[i]-floor(t[i]) == t[flag] - floor(t[flag]) && t[i] < t[flag]))
        {
          
          flag=i;
          
          
          
        }
    }
    norm_count[flag]++;
    total_norm++;
    t[flag]=0;
  }   
  } //correcting for rounding errors
else
{

  norm_count[max]=norm_count[max]-(total_norm-tableSize);
}

///////////////////////////////////////////////////////////

//printf("\n M=%ld N=%ld", total_norm,total);

//************building the table****************

//symbol start positions 
long *cumul;
alloc_vector_int(&cumul,maxValue+2);
//cumul[0];
for(i=1;i<maxValue;i++)
  cumul[i]=cumul[i-1]+norm_count[i-1];
cumul[maxValue]=tableSize+1;



  
   // used to define the step size to distribute symbols given the table size. 
    //Have some idea now. As symbols should be spread as far as possible from each other
    //ie (0,1/2,1/4,3/4,....), any sequence with step size n/2+n/4+n/8+...+odd_number
    // will approimate this distribution. Also, as the step size and n value are coprime
    // all states will be visited exactly once.

 //#define FSE_TABLESTEP(tableSize) ((tableSize>>1) + (tableSize>>3) + 3) in original code
//substitute 11 and get step=1283
//now, we distribute the symbols
long step, pos=0;
step=(tableSize>>1)+(tableSize>>3)+3;
long *tableSymbol;
alloc_vector_int(&tableSymbol,tableSize);
for(i=0;i<maxValue;i++) //iterating over symbols
{
  for(j=0;j<norm_count[i];j++)  //iterating for the number of normalized occurences
  {
    tableSymbol[pos]=i;
    pos=(pos+step)&(tableSize-1);
  }
 }


 
 


// we build the table. The table is an arrya of struct.Each symbol has 3 fields
 //minBits, minState-> if the state value is higher than this, use minBits+1
 // dstate gives the change in state after encoding

 // Build table 
unsigned tableU16[4200];
 
    {   long u; for (u=0; u<tableSize; u++) {
        long s = tableSymbol[u];   
        tableU16[cumul[s]++] = (unsigned) (tableSize+u);   
    }   }

struct symbolTransform symbolTT[maxValue];
total=0;



for(i=0;i<maxValue;i++)
{
  switch(norm_count[i])
  {
    case 0:
    symbolTT[i].dstate=0;
    break;
    case 1:
    symbolTT[i].dBits=(tableLog<<16) - (1<< (tableLog-1));

    symbolTT[i].dstate=total-1;
    total++;
    break;

    default: ;
    unsigned maxBits = (long)(tableLog - (long)mylog2(norm_count[i]-1));
        unsigned minStatePlus = (long)(norm_count[i]<< (maxBits));
        symbolTT[i].dstate = total - norm_count[i];
        symbolTT[i].dBits= (maxBits << 16) - minStatePlus;
        total +=  norm_count[i];
        
        
   

  }
  }


struct symbolTransform symbolT;
unsigned stateValue=tableSize;
unsigned nBitsOut;

unsigned out_buff[srcSize];





write_long_bitwise_2(tableLog, 4 ,buffer);
write_long_bitwise_2(num_sym, 8 ,buffer);
write_long_bitwise_2(srcSize, sizeof(long)*8 ,buffer);

//printf("%ld %ld \n",tableLog,num_sym );
if(num_sym*(tableLog+8) < 256 * tableLog)

{
  for(i=0;i<maxValue;i++)  
    if(norm_count[i])
      {
         write_long_bitwise_2(i, 8 ,buffer);
        write_long_bitwise_2(norm_count[i], tableLog ,buffer);
      }
  

}
else
  for(i=0;i<maxValue;i++)
    {
      write_long_bitwise_2(norm_count[i], tableLog ,buffer);
      //printf("%ld \n",norm_count[i] );
    }



for(i=srcSize-1;i>=0;i--)
{
  symbolT=symbolTT[in_buff[i]];
  //nBitsOut=stateValue<symbolT.minStatePlus?symbolT.minBits:symbolT.minBits+1;
  
    out_buff[i]=stateValue;
  nBitsOut = (unsigned)(stateValue+symbolT.dBits) >> 16;
  
  //encode state value push bits(state value)

  
  stateValue= tableU16[(long)((stateValue>>(nBitsOut))+symbolT.dstate)];
  
}
//out_buff[nx*ny-1]=stateValue;
//*finalState=stateValue;
write_long_bitwise_2(stateValue, tableLog ,buffer);

for(i=0;i<srcSize;i++)
  {
    
    symbolT=symbolTT[in_buff[i]];
    nBitsOut = (unsigned)(out_buff[i]+symbolT.dBits) >> 16;
    write_long_bitwise_2(out_buff[i],nBitsOut,buffer);

    
  }
}
/*--------------------------------------------------------------------------*/
void compress_image(char* input_filename, FILE* output_file, long q,long* nx_r,long* ny_r) {

long *qmap_small = 0;
long *qmap = 0;


long** image=0;
long** quantised=0;


  

long nx, ny;
long i,j;
long* in_buff; 

BITBUFFER buffer;
alloc_bitbuffer(&buffer,nx*ny*32);

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
  //printf("%ld %ld \n",nx,ny );
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
  clock_t ti=clock();
  //encode_adaptive_wnc(in_buff,nx*ny,q,0.3,16384,0,output_file);
  /*FSE_compress(in_buff,srcSize/4,&buffer);
   FSE_compress(in_buff+(srcSize/4),srcSize/4,&buffer);
    FSE_compress(in_buff+(srcSize/2),srcSize/4,&buffer);
     FSE_compress(in_buff+(3*srcSize/4),srcSize/4,&buffer);*/

  long blockSize= 65536;
  while(blockSize < srcSize/2)
  {
    blockSize= blockSize<<1;
  }

  //blockSize=blockSize>>1;
  long numBlocks= (long)(srcSize/blockSize);
  long rest= (long)(srcSize%blockSize);

  write_long_bitwise_2(blockSize,sizeof(long)*8,&buffer);
   write_long_bitwise_2(numBlocks,sizeof(long)*8,&buffer);
    write_long_bitwise_2(rest,sizeof(long)*8,&buffer);

    i=0;
     for(i=0;i<numBlocks;i++)
     {
      FSE_compress(in_buff+(i*blockSize),blockSize,&buffer);

     }

     if(rest)
      FSE_compress(in_buff+(i*blockSize),rest,&buffer);

 	//printf("here %ld \n", in_buff[k-1]);
//finding optimal value of tableLog
 	//minTableLog should be chosen from minBitsSrc and minBitsSymbols, as the source size > no.of symbols
 	// skip the computation and chose minBitsSrc as minBits, but 11>log(255)+2, and maxBitsSrc>11, we choose 11 as default tableLog
 	// by making some assumptions about srcSize and no. of symbols, the default value in the original code is 11, and here we skip
 	// the computations and just choose 11  (refer fse_compress -> optimalTableLog)



  
  bitbuffer_writefile(&buffer,output_file);
  ti=clock()-ti;
  double cpu_time=((double)ti)/CLOCKS_PER_SEC;
  printf("\n time to compress in ms: %f", cpu_time*1000);


// start encoding from the last symbol and go to the first symbol , but first we have to
//initialise state


  
}

void fse_decompress(long* out_buff, BITBUFFER *buffer)
{
  long tableLog,num_sym,srcSize;

  tableLog=bitbuffer_getbits(buffer,4);
  num_sym=bitbuffer_getbits(buffer,8);
  srcSize=bitbuffer_getbits(buffer,sizeof(long)*8);

  //printf("%ld %ld \n",tableLog,num_sym );

  long tableSize=1<<tableLog;
  
  long *norm_count;
  alloc_vector_int(&norm_count,256);
  struct decode_symbol
  {
    unsigned newState;
    unsigned symbol;
    unsigned nBits;
  } decode_table[tableSize];

  long pos=0;
  long step=(tableSize>>1)+(tableSize>>3)+3;
  long *tableSymbol;
  alloc_vector_int(&tableSymbol,tableSize);
  unsigned i,j;

  

  

  for(i=0;i<256;i++)
    norm_count[i]=0;
  
  if(num_sym*(tableLog+8) < 256 * tableLog)
  {
    for(i=0;i<num_sym;i++)  
      {
        long k= bitbuffer_getbits(buffer,8);
        norm_count[k]=bitbuffer_getbits(buffer,tableLog);
        for(j=0;j<norm_count[k];j++)  //iterating for the number of normalized occurences
    {
      decode_table[pos].symbol=k;
      pos=(pos+step)&(tableSize-1);
    }
      }
  }
  else
  {


  for(i=0;i<256;i++)
    {
      norm_count[i]=bitbuffer_getbits(buffer,tableLog);
      //printf("%ld \n", norm_count[i]);
      for(j=0;j<norm_count[i];j++)  //iterating for the number of normalized occurences
    {
      decode_table[pos].symbol=i;
      pos=(pos+step)&(tableSize-1);
    }
      //printf("%ld -> %ld \n", i, norm_count[i]);
    }
  }

  

  //build decoding table

  for(i=0;i<tableSize;i++)
  {
    unsigned sym=decode_table[i].symbol;
    unsigned nextState= norm_count[sym]++;
    decode_table[i].nBits = (unsigned)(tableLog - mylog2(nextState));
    //if(decode_table[i].nBits>tableLog)
      //decode_table[i].nBits=tableLog;
    
    decode_table[i].newState=(unsigned)((nextState<< decode_table[i].nBits) - tableSize);
    if(decode_table[i].newState==tableSize)
      {decode_table[i].newState=(unsigned)(((nextState-1)<< decode_table[i].nBits) - tableSize); }
    //printf("state->%ld symbol->%ld nBits->%ld nextState->%ld newState->%ld \n", i, decode_table[i].symbol, decode_table[i].nBits, nextState, decode_table[i].newState);
    
  }

  //unsigned initialState = 


  //printf("%ld \n",initialState );

  unsigned stateValue=bitbuffer_getbits(buffer, tableLog);;
  //printf("%ld \n", stateValue);
  struct decode_symbol stateInfo;
  long nextBits;
  long k=0;
    
  while(k<srcSize)
  {
      
      stateInfo=decode_table[stateValue];
    out_buff[k]=stateInfo.symbol;

    //if(k>=nx*ny-500)
    //printf("%ld->%ld %ld %ld",k,out_buff[k], stateValue,stateInfo.nBits );
    
    //nextBits=read_long_bitwise_2(&buffer,stateInfo.nBits);
    nextBits=bitbuffer_getbits(buffer,stateInfo.nBits);
  
    stateValue=(stateInfo.newState + nextBits) ;
     k++;

  }
}

/*--------------------------------------------------------------------------*/
void decompress_image(FILE* input_file, char* output_filename, long nx, long ny)
{
	

	// first step in decoding is building the decoding table. For now, we assume we have the norm_count values
	// but later we have to write it in the file and read from it.

	//first spread the symbols


  //long out_buff[nx*ny];
  //decode_adaptive_wnc(input_file,nx*ny,256,0.3,16384,0,out_buff);
  long i,j;

  BITBUFFER buffer;
  alloc_bitbuffer(&buffer,(nx*ny)<<5);
  
  bitbuffer_loadfile(&buffer, input_file);
  long out_buff[nx*ny];
  /*long out_buff1[nx*ny/4];
  long out_buff2[nx*ny/4];
  long out_buff3[nx*ny/4];
  long out_buff4[nx*ny/4];*/
   clock_t t,t_tot=0;
    t=clock();
	 
   /*fse_decompress(out_buff,&buffer);
   fse_decompress(out_buff + (nx*ny/4),&buffer);
   fse_decompress(out_buff + (nx*ny/2),&buffer);
   fse_decompress(out_buff + (3*nx*ny/4),&buffer);*/

    long blockSize=bitbuffer_getbits(&buffer,sizeof(long)*8);
   long numBlocks=bitbuffer_getbits(&buffer,sizeof(long)*8);
    long rest=bitbuffer_getbits(&buffer,sizeof(long)*8);

    printf("%ld %ld %ld \n",blockSize,numBlocks,rest );

    i=0;
     for(i=0;i<numBlocks;i++)
     {
      fse_decompress(out_buff+(i*blockSize),&buffer);
     }

     if(rest)
      fse_decompress(out_buff+(i*blockSize),&buffer);


  t_tot+=clock()-t;
  double cpu_time=((double)t_tot)/CLOCKS_PER_SEC;
  printf("\n time to decompress: %f ms", cpu_time*1000);

  /*for(i=0;i<nx*ny/4;i++)
    {
      out_buff[i]=out_buff1[i];
      out_buff[(nx*ny/4)+i]=out_buff2[i];
      out_buff[(nx*ny/2)+i]=out_buff3[i];
      out_buff[(3*nx*ny/4)+i]=out_buff4[i];
    }*/



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

  FILE* output_file=0;

  FILE* input_file=0;
long nx,ny;

//double cpu_time;


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
   

  //unsigned norm_count[256];
  //long tableSymbol[2048];
  //t=clock();
  output_file = fopen(compressed_filename, "w");
  compress_image (input_filename,output_file,q,&nx,&ny);
  //write_long_bitwise(2693,4,output_file);
  //fclose(input_file);
  fclose(output_file);
  //t=clock()-t;
  //double cpu_time=((double)t)/CLOCKS_PER_SEC;
  //printf("\n time to compress: %f", cpu_time);

  
  input_file = fopen(compressed_filename, "r");
  decompress_image(input_file,output_filename,nx,ny);
  fclose(input_file);
  

  
  
  return 0;
}
