#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <string.h>

#include "tv.h"


/*! \file main.c
  
  \brief This file contains methods for executing the program
 */

/*!
  This function executes the different parts of the program according
  to the input options.
 */

int main (int argc, char **argv)
{
  
  if (argc < 2) {
    fprintf (stderr, "\n\n");
    fprintf (stderr, "Usage\n");
    fprintf (stderr, " 0: \n");

    fprintf (stderr, "\n\n");

    exit (0);
  }

  switch (atoi(argv[1])) {
    
    /*! 0: 
      
     */
  case 0 : {
    


    break;
  }


  default : {

    

    break;
  }
  }


  return 0;
}

/*! \mainpage TV DOCUMENTATION

  \author 
  Carlos Vera-Ciro \n
  Astronomy Department \n
  Madison University \n
  475 N Charter St  Madison, WI 53706 \n
  ciro@astro.wisc.edu \n\n


  \version
  0.alpha
  
*/

/*! \page Makefile Makefile
  
  - <b> SIMULATION=<sim> </b>: Selects the simulation that is going to
  be analyzed.
  

*/

