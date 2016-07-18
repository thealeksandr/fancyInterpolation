//
//  akima.c
//  spline
//
//  Created by Александр Никифоров on 7/18/16.
//  Copyright © 2016 thealeksandr. All rights reserved.
//

#include <stdio.h>

/***************************************************************************
 * aspline.c -- aspline does an Akima spline interpolation.                *
 *                                                                         *
 * spline is (c) David Frey, 1996, 1998, 1999, 2002, 2009                  *
 *                                                                         *
 * This program is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU General Public License as published by the   *
 * Free Software Foundation; either version 2 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * This program is distributed in the hope that it will be useful, but     *
 * WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * General Public License for more details.                                *
 *                                                                         *
 * You should have received a copy of the GNU General Public License       *
 * along with this program; if not, write to the Free Software Foundation, *
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA            *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include <errno.h>

#define DXSTEP          1.0
#define INTERVALLS      100

const char *progname;

double *x, *y;           /* the values to interpolate */

int warranty(void);
int usage(void);
void insertpoint(const char *inname, double xv, double yv, int *n);
void computelimits(int s, int e, double *llimit, double *ulimit);
void calcspline(int n, int d, double xstep, double llimit, double ulimit,
                int calclimit);
void readandcalcspline(FILE *infile, const char *inname,
                       int d, double xstep, double llimit, double ulimit,
                       int calclimit, int verbose, int *s);

static struct option const long_options[] = {
    {"auto",     optional_argument, 0, 'a'},
    {"help",     no_argument,       0, 'h'},
    {"help",     no_argument,       0, '?'},
    {"llimit",   required_argument, 0, 'l'},
    {"ulimit",   required_argument, 0, 'u'},
    {"points",   required_argument, 0, 'n'},
    {"verbose",  no_argument,       0, 'v'},
    {"version",  no_argument,       0, 'V'},
    {"license",  no_argument,       0, 'W'},
    {"warranty", no_argument,       0, 'W'},
    {0,          0,                 0, 0  }
};

int warranty(void)
/* say that there is no warranty */
{
    
    printf("This program is free software; you can redistribute it and/or "\
           "modify it\n"\
           "under the terms of the GNU General Public License as published "\
           "by the\n"\
           "Free Software Foundation; either version 2 of the License, or "\
           "(at your\n"\
           "option) any later version.\n"\
           "\n"\
           "This program is distributed in the hope that it will be useful, "\
           "but\n"\
           "WITHOUT ANY WARRANTY; without even the implied warranty of "\
           "MERCHANTABILITY\n"\
           "or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public "\
           "License\n"\
           "for more details.\n"\
           "\n"\
           "You should have received a copy of the GNU General Public License "\
           "along\n"\
           "with this program; if not, write to the Free Software Foundation, "\
           "Inc.,\n"\
           "59 Temple Place Suite 330, Boston, MA 02111-1307 USA.\n");
    exit(0);
}

int usage(void)
/* print a short usage statement. */
{
    fprintf(stderr,
            "usage: %s [-hv][-a [xstep]][-l llimit][-u ulimit][-n points] "\
            "{data-file(s)}\n\n",progname);
    fprintf(stderr,
            "  -a, --auto      supply abscissas automatically.\n"\
            "  -l, --llimit    set lower x-limit.\n"\
            "  -u, --ulimit    set upper x-limit.\n"\
            "  -n, --points    space output points so that approx. "\
            "n intervals\n"\
            "                  occur between the lower and upper limit.\n"\
            "  -h, --help      display this help and exit.\n"\
            "  -v, --verbose   be verbose.\n"\
            "  -V, --version   display version and copyright information and "\
            "exit.\n"\
            "  -W, --warranty  display licensing terms.\n");
    
    exit(0);
}

void insertpoint(const char *inname, double xv, double yv, int *n)
/* insert or drop the point if duplicate */
{
    int k;
    
    /* search the point */
    k=0;
    while ((k < *n) && (xv > x[k])) k++;
    
    if (k == *n) { /* new point */
        /* We malloc the first point explicitly, on some architectures an
         * malloc(NULL,...) != realloc(...)                               */
        if (*n == 0) {
            x=xmalloc(sizeof(double));
            y=xmalloc(sizeof(double));
        } else {
            x=xrealloc(x, (*n+1)*sizeof(double));
            y=xrealloc(y, (*n+1)*sizeof(double));
        }
        
        x[*n]=xv; y[*n]=yv; (*n)++;
    } else if (fpclassify(xv - x[k]) == FP_ZERO) { /* duplicate */
        fprintf(stderr,
                "%s:%s: duplicate abscissas found, sample (%g,%g), dropped.\n",
                progname, inname, xv,yv);
    } else { /* intermediate point. Insert it */
        x=xrealloc(x, (*n+1)*sizeof(double)); y=xrealloc(y, (*n+1)*sizeof(double));
        
        /* Move the points above */
        memmove(&x[k+1],&x[k], (size_t)((*n-k)*sizeof(double)));
        memmove(&y[k+1],&y[k], (size_t)((*n-k)*sizeof(double)));
        
        x[k]=xv; y[k]=yv; (*n)++;
    }
}

void computelimits(int s, int e, double *llimit, double *ulimit)
/* Compute minimum and maximum of the x-values */
{
    int i;
    
    *llimit=FLT_MAX; *ulimit=-FLT_MAX;
    for (i=s; i < e; i++) {
        if (x[i] < *llimit) *llimit=x[i];
        if (x[i] > *ulimit) *ulimit=x[i];
    }
}

void calcspline(int n, int d, double xstep, double llimit, double ulimit,
                int calclimit)
/* Calculate the spline and output the result on stdout */
{
    if (n > 0) { /* we have data to process */
        if (n == 1) { /* but only one point :( */
            printf("%g %g\n", x[0], y[0]);
        } else if (n == 2) { /* wow, two points! */
            double dx, dy, m;
            int i;
            
            dx=x[1]-x[0]; dy=y[1]-y[0];
            
            /* 3rd step: output the coefficients for the subintervalls i=2..n-4 */
            /* calculate the intermediate values */
            
            if (calclimit) computelimits(0,n,&llimit, &ulimit);
            xstep=(ulimit-llimit)/d;
            
            m=dy/dx; /* dx != 0.0 asserted by insertpoint */
            for (i=0; i <= d; i++)
                printf("%g %g\n", x[0]+i*xstep, y[0]+i*m*xstep);
        } else { /* now we have something to compute */
            int i, p;                  /* n=no of points */
            double *dx, *dy;           /* x[i+1]-x[i], resp. y[i+1]-y[i] */
            double *m, *t;             /* the slopes */
            double *C, *D;             /* the coefficients */
            
            /* Leading extrapolation points,
             actual values will be filled in later */
            
            x=(double *)xrealloc(x, (n+4) * sizeof(double));
            y=(double *)xrealloc(y, (n+4) * sizeof(double));
            
            /* Move x[0] to x[2] */
            memmove(&x[2],&x[0], (size_t)(n*sizeof(double)));
            memmove(&y[2],&y[0], (size_t)(n*sizeof(double)));
            n += 4;
            
            /* calculate coefficients of the spline (Akima interpolation itself) */
            
            /* allocate the arrays
             * NB: There are some unused locations (dx[0]),
             *     but this faciliates the indexing.
             */
            
            dx=(double *)xmalloc(n * sizeof(double));
            dy=(double *)xmalloc(n * sizeof(double));
            m=(double *)xmalloc(n * sizeof(double));
            t=(double *)xmalloc(n * sizeof(double));
            
            /* a) Calculate the differences and the slopes m[i]. */
            
            for (i=2; i < n-3; i++) {
                dx[i]=x[i+1] - x[i]; dy[i]=y[i+1] - y[i];
                m[i] =dy[i]/dx[i]; /* dx != 0, asserted by insertpoint() */
            }
            
            /* b) interpolate the missing points: */
            
            x[1]=x[2] + x[3] - x[4];           dx[1]=x[2] - x[1];
            y[1]=dx[1]*(m[3] - 2*m[2]) + y[2]; dy[1]=y[2] - y[1];
            m[1]=dy[1]/dx[1];
            
            x[0]=2*x[2] - x[4];                dx[0]=x[1] - x[0];
            y[0]=dx[0]*(m[2] - 2*m[1]) + y[1]; dy[0]=y[1] - y[0];
            m[0]=dy[0]/dx[0];
            
            x[n-2]=x[n-3] + x[n-4] - x[n-5];
            y[n-2]=(2*m[n-4] - m[n-5])*(x[n-2] - x[n-3l]) + y[n-3];
            
            x[n-1]=2*x[n-3] - x[n-5];
            y[n-1]=(2*m[n-3] - m[n-4])*(x[n-1] - x[n-2]) + y[n-2];
            
            for (i=n-3; i < n-1; i++) {
                dx[i]=x[i+1] - x[i]; dy[i]=y[i+1] - y[i];
                m[i]=dy[i]/dx[i];
            }
            
            /* the first x slopes and the last y ones are extrapolated: */
            
            t[0]=0.0; t[1]=0.0;  /* not relevant */
            for (i=2; i < n-2; i++) {
                double num, den;
                
                num=fabs(m[i+1] - m[i])*m[i-1] + fabs(m[i-1] - m[i-2])*m[i];
                den=fabs(m[i+1] - m[i]) + fabs(m[i-1] - m[i-2]);
                
                if (fpclassify(den) != FP_ZERO) t[i]=num/den;
                else                            t[i]=0.0;
            }
            
            /* c) Allocate the polynom coefficients */
            
            /* A=y_i
             B=t_i */
            C=(double *)xmalloc(n * sizeof(double));
            D=(double *)xmalloc(n * sizeof(double));
            
            for (i=2; i < n-2; i++) {
                /* A[i]=y[i];
                 B[i]=t[i]; */
                C[i]=(3*m[i] - 2*t[i] - t[i+1])/dx[i];
                D[i]=(t[i] + t[i+1] - 2*m[i])/(dx[i]*dx[i]);
            }
            
            /* 3rd step: output the coefficients for the subintervalls i=2..n-4 */
            /* calculate the intermediate values */
            
            if (calclimit) computelimits(2,n-2,&llimit, &ulimit);
            xstep=(ulimit-llimit)/d;
            
            p=2;
            double xv;
            for(xv=llimit; xv<ulimit+xstep; xv += xstep) {
                while (xv >= x[p]) {
                    printf("%g %g\n", x[p], y[p]); p++;
                }
                
                
                /* skip the next interpolated point if it's too close to the current
                 point */
                if (((xv - x[p-1]) > xstep/100.0) &&
                    ((x[p] - xv) > xstep/100.0)) {
                    double xd=(xv-x[p-1]);
                    printf("%g %g\n", xv,
                           y[p-1] + (t[p-1] + (C[p-1] + D[p-1]*xd)*xd)*xd);
                }
            }
            
            free(dx); free(dy); free(m); free(t); free(C); free(D);
        }
    }
}

void readandcalcspline(FILE *infile, const char *inname,
                       int d, double xstep, double llimit, double ulimit,
                       int calclimit, int verbose, int *s)
/* s is the dataset number */
{
    int n, l;
    char *line;
    
    if (infile != stdin) {
        infile=fopen(inname, "r");
        if (infile == NULL) {
            fprintf(stderr,
                    "%s: can't open \"%s\": %s.\n",
                    progname, inname, strerror(errno));
            return;
        }
    }
    
    if (verbose) printf("# %s\n", inname);
    
    /* Read the points and build an array with the samples
     (using insertpoint()). */
    
    n=0; l=0; x=NULL; y=NULL; line=NULL;
    do
    {
        char *il;
        
        if (line != NULL) free(line);
        
        line=gettextline(infile);
        
        il=line; l++;
        if (il != NULL)
            while (isspace(*il)) il++;
        
        if ((line == NULL) || (il[0] == '\0')) {
            if (n > 0) {
                if (*s > 0) printf("\n");
                calcspline(n, d, xstep, llimit, ulimit, calclimit); (*s)++;
                free(x); free(y); x=y=NULL; n=0;
            }
        } else if (il[0] != '#') {
            int k;                   /* number of values */
            double f1, f2;            /* doubleing point numbers */
            
            /* read point */
            k=sscanf(line, "%lf %lf", &f1, &f2);
            
            if (k < 1)
                fprintf(stderr, "%s:%s:%d: parse error at '%s'\n",
                        progname, inname, l, line);
            else if (k==1)
                insertpoint(inname, llimit+n*xstep,f1,&n); /* Auto xscale */
            else
                insertpoint(inname,f1,f2,&n);
        }
    }
    while (line != NULL);
    
    if (infile != stdin) fclose(infile);
}

int main(int argc, char *argv[])
{
    int c;
    int i,d;
    double xstep, llimit, ulimit;  /* x-spacing, lower limit, upper limit */
    int calclimit, verbose;       /* calculate limits, be verbose */
    
    progname=strrchr(argv[0], '/');
    if (progname == NULL) progname=argv[0];
    else                  progname++;
    
    /*
     * These are mostly preliminary assignments, the actual values will
     * be calculated when evaluating the spline expression depending on
     */
    
    verbose=0; calclimit=1; xstep=DXSTEP; llimit=ulimit=0.0; d=INTERVALLS;
    while ((c = getopt_long(argc, argv, "a:?hl:n:u:vVW",
                            long_options, (int *)0)) != EOF) {
        switch (c) {
            case 0:   break;
            case 'a': xstep=strtod(optarg, NULL);
                if (fpclassify(xstep) == FP_ZERO) xstep=DXSTEP;
                break;
            case 'l': llimit=strtod(optarg, NULL); calclimit=0; break;
            case 'h':
            case '?': usage(); break;
            case 'n': d=atoi(optarg); break;
            case 'u': ulimit=strtod(optarg, NULL); calclimit=0; break;
            case 'v': verbose++; break;
            case 'V': fprintf(stderr, "%s version %s\n"
                              "Akima spline interpolation (c) 1996-1998 David Frey.\n"
                              "This is free software; see the GNU General Public Licence version 2 or "
                              "later\n"
                              "for copying conditions.\n"
                              "There is NO warranty.  See %s --licence for details.\n",
                              progname, VERSION, progname);
                return 0; break;
            case 'W': warranty(); break;
            default : break;
        }
    }
    
    if (ulimit < llimit) {
        double temp;
        
        temp=ulimit; ulimit=llimit; llimit=temp;
    }
    
    if (argc > optind) {
        int s;
        
        s=0;
        for (i=optind; i < argc; i++) {
            readandcalcspline(NULL, argv[i],
                              d, xstep, llimit, ulimit, calclimit, verbose,
                              &s);
            
            if ((s > 0) && (i < argc-1)) printf("\n");
            s++;
        }
    } else {
        int s;
        
        s=0;
        readandcalcspline(stdin, "stdin",
                          d, xstep, llimit, ulimit, calclimit, verbose,
                          &s);
    }
    return 0;
}
