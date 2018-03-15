#include "renderer.h"
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <stdio.h>

#include <png.h>
#include <stdlib.h>
#include <stdint.h>


void generateRaysFromCamera(camera cam, tracingRay * rays, unsigned int numRaysRows, unsigned int numRaysCols)
{
   printf("starting routine\n");
   struct renderman_S renderman;
   renderman.camU = gsl_vector_alloc(3);
   renderman.camV = gsl_vector_alloc(3);
   renderman.camW = gsl_vector_alloc(3);
   renderman.rot = gsl_matrix_alloc(4,4);
   renderman.trans = gsl_matrix_alloc(4,4);
   renderman.boxpoints = gsl_matrix_alloc(4,8);
   renderman.renderBoxpoints = gsl_matrix_alloc(4,8);
   renderman.coordinateTransform = gsl_matrix_alloc(4,4);
   gsl_matrix * inverseTransform = gsl_matrix_alloc(4,4);
   gsl_matrix * originPoint = gsl_matrix_alloc(4,1);
   gsl_matrix * unitDirectionVector = gsl_matrix_alloc(4,1);
   gsl_matrix * result = gsl_matrix_alloc(4,1);
   
   //have camera pose and euler angles, need resolution and size of pixels 
   //and have distance of film from pinhole as zoom

   //parameterize length and height on indices seperated in space by ray density distance
   //let raydistance be .01
   double raydistance = .01;
   double centerx = raydistance * ((double)numRaysRows)/2.0;
   double centery = raydistance * ((double)numRaysCols)/2.0;
   //start at {(centerx - raydistance*numRays/2) , (centery - raydistance*numRays/2)}
   //loop by raydistance, nested dims

   //transform from camera imaging plane to world coord

  //rotate by pan and tilt(rotation around Z/W)
  gsl_vector_set_basis(renderman.camU, 0);
  gsl_vector_set_basis(renderman.camV, 1);
  gsl_vector_set_basis(renderman.camW, 2);
  gsl_blas_drot(renderman.camU, renderman.camV, cos(cam.roll), sin(cam.roll));
  gsl_blas_drot(renderman.camU, renderman.camW, cos(cam.pan), sin(cam.pan));
  gsl_blas_drot(renderman.camV, renderman.camW, cos(cam.tilt), sin(cam.tilt));
  
  printf("forming transform matrices\n");
  //form rotation and translation matrices
  for (unsigned int i=0; i<3; i++)
  {
    gsl_matrix_set(renderman.trans, i, 3, gsl_vector_get(cam.pose, i));
    gsl_matrix_set(renderman.rot,   i, 0, gsl_vector_get(renderman.camU, i));
    gsl_matrix_set(renderman.rot,   i, 1, gsl_vector_get(renderman.camV, i));
    gsl_matrix_set(renderman.rot,   i, 2, gsl_vector_get(renderman.camW, i));
  }
  gsl_matrix_set(renderman.rot, 3,3,1);
  for (unsigned int i=0; i<4; i++)
       {gsl_matrix_set(renderman.trans, i,i,1);}

  printf("\nrotation matrix: \n");
  for (unsigned int i=0; i<4; i++)
  printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.rot, i,0),
                              gsl_matrix_get(renderman.rot, i,1),
                              gsl_matrix_get(renderman.rot, i,2), 
                              gsl_matrix_get(renderman.rot, i,3));

  printf("\ntranslation matrix: \n");
  for (unsigned int i=0; i<4; i++)
  printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.trans, i,0),
                              gsl_matrix_get(renderman.trans, i,1),
                              gsl_matrix_get(renderman.trans, i,2), 
                              gsl_matrix_get(renderman.trans, i,3));

  //calculate total transformation matrix
  printf("calculating total tf\n");
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,  1, renderman.trans, renderman.rot,
                 0, renderman.coordinateTransform);

  printf("\ntransform matrix: \n");
  for (unsigned int i=0; i<4; i++)
  printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.coordinateTransform, i,0), 
                              gsl_matrix_get(renderman.coordinateTransform, i,1),
                              gsl_matrix_get(renderman.coordinateTransform, i,2), 
                              gsl_matrix_get(renderman.coordinateTransform, i,3));

  //calculate inverse of this transorm
  printf("calculating inverse\n");
  gsl_permutation * myPerm =  gsl_permutation_alloc (4);
  int signum;
  printf("calculating inverse2\n");
  gsl_linalg_LU_decomp(renderman.coordinateTransform, myPerm, &signum);
  printf("calculating inverse3\n");
  gsl_linalg_LU_invert(renderman.coordinateTransform, myPerm, inverseTransform);
  gsl_permutation_free(myPerm);
  //now have start point in each loop, define each ray 
  //as starting there and having a unit vector pointing towards the point zoom distance away from
  //the centerx, centery point and normal to the imaging plane
  // (so vector from point to centerx,y + normal vector. Then normalized)
  printf("going to generate rays\n");
  //calculate normal point (tip of normal vector) of imaging plane with zoom length
  // (thats 0,0,0)
  double normalx = 0;
  double normaly = 0;
  double normalz = 0;

  for (int row=0; row<numRaysRows; row++)
  {
     for (int col=0; col<numRaysCols; col++)
     {
        //find point on imaging plane in camera coord
        double filmy = row * raydistance - centery;
        double filmx = col * raydistance - centerx;
        double filmz = -cam.zoom;
        //take (normal point - film point) and normalize to get unit direction 
        // vector in camera coord
        double magnitude = sqrt(filmy * filmy + filmx * filmx + cam.zoom*cam.zoom);
        
        double udvU = -filmy/magnitude;
        double udvV = -filmx/magnitude;
        double udvW = cam.zoom/magnitude;
        //apply inverse tf to film point to get origin in world coord
        gsl_matrix_set(originPoint, 0,0, filmx);
        gsl_matrix_set(originPoint, 1,0, filmy);
        gsl_matrix_set(originPoint, 2,0, filmz);
        gsl_matrix_set(originPoint, 3,0, 1);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, inverseTransform, 
                       originPoint, 0, result);
        
        rays[row + numRaysRows * col].origin[0] = gsl_matrix_get(result, 0,0);
        rays[row + numRaysRows * col].origin[1] = gsl_matrix_get(result, 1,0);
        rays[row + numRaysRows * col].origin[2] = gsl_matrix_get(result, 2,0);
        //apply inverse tf to unit direction vector to get udv in world coord
        gsl_matrix_set(unitDirectionVector, 0,0, udvU);
        gsl_matrix_set(unitDirectionVector, 1,0, udvV);
        gsl_matrix_set(unitDirectionVector, 2,0, udvW);
        gsl_matrix_set(unitDirectionVector, 3,0, 1);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, renderman.rot, 
                       unitDirectionVector, 0, result);

        rays[row + numRaysRows * col].direction[0] = gsl_matrix_get(result, 0,0);// - rays[row + numRaysRows * col].origin[0];
        rays[row + numRaysRows * col].direction[1] = gsl_matrix_get(result, 1,0);// - rays[row + numRaysRows * col].origin[1];
        rays[row + numRaysRows * col].direction[2] = gsl_matrix_get(result, 2,0);// - rays[row + numRaysRows * col].origin[2];
        //rays[row + numRaysRows * col].direction[0] = udvU;
        //rays[row + numRaysRows * col].direction[1] = udvV;
        //rays[row + numRaysRows * col].direction[2] = udvW;
     }
  }
  
  printf("freeing allocated memory");
  gsl_vector_free (renderman.camU);
  gsl_vector_free (renderman.camV);
  gsl_vector_free (renderman.camW);
  gsl_matrix_free (renderman.rot); 
  gsl_matrix_free (renderman.trans);
  gsl_matrix_free (renderman.boxpoints);
  gsl_matrix_free (renderman.renderBoxpoints);
  gsl_matrix_free (inverseTransform);
  gsl_matrix_free (originPoint);
  gsl_matrix_free (unitDirectionVector);
  gsl_matrix_free (result);
  
}

//Following routines borrowed from https://www.lemoda.net/c/write-png/
/* A coloured pixel. */

typedef struct {
    uint8_t red;
    uint8_t green;
    uint8_t blue;
} pixel_t;

/* A picture. */
    
typedef struct  {
    pixel_t *pixels;
    size_t width;
    size_t height;
} bitmap_t;
    
/* Given "bitmap", this returns the pixel of bitmap at the point 
   ("x", "y"). */

static pixel_t * pixel_at (bitmap_t * bitmap, int x, int y)
{
    return bitmap->pixels + bitmap->width * y + x;
}
    
/* Write "bitmap" to a PNG file specified by "path"; returns 0 on
   success, non-zero on error. */

static int save_png_to_file (bitmap_t *bitmap, const char *path)
{
    FILE * fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    size_t x, y;
    png_byte ** row_pointers = NULL;
    /* "status" contains the return value of this function. At first
       it is set to a value which means 'failure'. When the routine
       has finished its work, it is set to a value which means
       'success'. */
    int status = -1;
    /* The following number is set by trial and error only. I cannot
       see where it it is documented in the libpng manual.
    */
    int pixel_size = 3;
    int depth = 8;
    
    fp = fopen (path, "wb");
    if (! fp) {
        goto fopen_failed;
    }

    png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        goto png_create_write_struct_failed;
    }
    
    info_ptr = png_create_info_struct (png_ptr);
    if (info_ptr == NULL) {
        goto png_create_info_struct_failed;
    }
    
    /* Set up error handling. */

    if (setjmp (png_jmpbuf (png_ptr))) {
        goto png_failure;
    }
    
    /* Set image attributes. */

    png_set_IHDR (png_ptr,
                  info_ptr,
                  bitmap->width,
                  bitmap->height,
                  depth,
                  PNG_COLOR_TYPE_RGB,
                  PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_DEFAULT,
                  PNG_FILTER_TYPE_DEFAULT);
    
    /* Initialize rows of PNG. */

    row_pointers = png_malloc (png_ptr, bitmap->height * sizeof (png_byte *));
    for (y = 0; y < bitmap->height; y++) {
        png_byte *row = 
            png_malloc (png_ptr, sizeof (uint8_t) * bitmap->width * pixel_size);
        row_pointers[y] = row;
        for (x = 0; x < bitmap->width; x++) {
            pixel_t * pixel = pixel_at (bitmap, x, y);
            *row++ = pixel->red;
            *row++ = pixel->green;
            *row++ = pixel->blue;
        }
    }
    
    /* Write the image data to "fp". */

    png_init_io (png_ptr, fp);
    png_set_rows (png_ptr, info_ptr, row_pointers);
    png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* The routine has successfully written the file, so we set
       "status" to a value which indicates success. */

    status = 0;
    
    for (y = 0; y < bitmap->height; y++) {
        png_free (png_ptr, row_pointers[y]);
    }
    png_free (png_ptr, row_pointers);
    
 png_failure:
 png_create_info_struct_failed:
    png_destroy_write_struct (&png_ptr, &info_ptr);
 png_create_write_struct_failed:
    fclose (fp);
 fopen_failed:
    return status;
}

/* Given "value" and "max", the maximum value which we expect "value"
   to take, this returns an integer between 0 and 255 proportional to
   "value" divided by "max". */

static int pix (int value, int max)
{
    if (value < 0) {
        return 0;
    }
    return (int) (256.0 *((double) (value)/(double) max));
}
//end of code from https://www.lemoda.net/c/write-png/

//following code taken from programmingalgorithms.com
struct RGB
{
	unsigned char R;
	unsigned char G;
	unsigned char B;
};

struct HSL
{
	int H;
	float S;
	float L;
};

float HueToRGB(float v1, float v2, float vH)
{
	if (vH < 0)
		vH += 1;

	if (vH > 1)
		vH -= 1;

	if ((6 * vH) < 1)
		return (v1 + (v2 - v1) * 6 * vH);

	if ((2 * vH) < 1)
		return v2;

	if ((3 * vH) < 2)
		return (v1 + (v2 - v1) * ((2.0f / 3) - vH) * 6);

	return v1;
}

struct RGB HSLToRGB(struct HSL hsl) {
	struct RGB rgb;

	if (hsl.S == 0)
	{
		rgb.R = rgb.G = rgb.B = (unsigned char)(hsl.L * 255);
	}
	else
	{
		float v1, v2;
		float hue = (float)hsl.H / 360;

		v2 = (hsl.L < 0.5) ? (hsl.L * (1 + hsl.S)) : ((hsl.L + hsl.S) - (hsl.L * hsl.S));
		v1 = 2 * hsl.L - v2;

		rgb.R = (unsigned char)(255 * HueToRGB(v1, v2, hue + (1.0f / 3)));
		rgb.G = (unsigned char)(255 * HueToRGB(v1, v2, hue));
		rgb.B = (unsigned char)(255 * HueToRGB(v1, v2, hue - (1.0f / 3)));
	}

	return rgb;
}
//end of code from programmingalgorithms.com

//render is ours again, it uses the above borrowed code
void render(camera cam, tracingRay * rays, unsigned int numRaysRows, unsigned int numRaysCols, const char* filename)
{
  bitmap_t image;
  int x;
  int y;

  image.width = numRaysRows;
  image.height = numRaysCols;
  //for now, each ray is a pixel
  struct HSL data;
  struct RGB values;
  image.pixels = calloc(image.width * image.height, sizeof(pixel_t));
  for (y=0; y <image.height; y++)
  {
     for (x=0; x<image.width; x++)
     {
       pixel_t * pix = pixel_at(&image, x, y);
       data.H =  rays[y + numRaysRows * x].hue;
       data.S = rays[y + numRaysRows * x].sat;
       data.L = rays[y + numRaysRows * x].lightness;
       values = HSLToRGB(data);
       pix->red = values.R;
       pix->green = values.G;
       pix->blue = values.B;
     } 
  }

  save_png_to_file(&image, filename);
  free (image.pixels);
}