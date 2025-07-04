/*
 *	StarGen Main Routine 
 *
 *	This file provides the main command-line interface to StarGen.
 *	Other platform-specific UIs can be created by duplicating its
 *	general functionality and then calling stargen(), whose API is
 *	defined in stargen.h
 *
 *	$Id: main.c,v 1.44.0006 2025/07/04 08.31 $ 
 */

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<time.h>
#include	<ctype.h>

#include	"Dumas.h"

#ifdef THINK_C
#define macintosh 1
#endif

#ifdef macintosh
#include	<console.h>
#include	<unix.h>
#else
#include	<sys/types.h>
#endif

#ifdef MSDOS
#include	<stddef.h>
#include	<stdlib.h>
#include	<float.h>
#endif

#ifdef WIN32
#pragma warning (disable: 4048) // Hush compiler - don't complain about storing specific arrays in generic pointers
#endif

#include	"const.h"
#include	"structs.h"
#include	"stargen.h"

//				extern double migratek;
// JN debug migration
double migratek=1.0;
// filter out too near planets, experimental buggy
int USE_FILTERING = 0;

double BHILL_CRITERION = 1.0; // eikä long double

int FILTER_ASTEROIDS = 0; // Set to 1 to remove asteroids, 0 to keep them
int USE_HILL = 1;         // Set to 1 to use Hill criterion, 0 for Orbital Period


// some basic params of cloud
extern long double dust_density_coeff; //2e-3
extern long double cloud_eccentricity; //0.2
extern long double nearestk; //0.3
extern long double farthestk; // 50
extern long double diskradiusk; // 200 
extern long double alphaa; // ALPHA 5?
extern long double gammaa; // N 3?
extern long double gasdust;
extern long double bee;

extern int use_exponent_disk; // use own disk surface density profile, alpha is exponent

extern int order_to_resonances;
extern long double base_resonance;

extern int use_own_luminosity;
extern long double par_luminosity;


/*
 *		StarGen supports private catalogs. The two here are ones I am using
 *		for debuggery. They may go away.
 */

#define EM(x)		(x)/SUN_MASS_IN_EARTH_MASSES
#define AVE(x,y)	((x+y)/2.)

				/*  No 	Orbit	Eccen. 	Tilt	Mass		Giant?	Dust Mass	Gas */
planets sphinx3   ={ 4,	3.0,	0.046,	10.5,	EM(2.35),	FALSE,	EM(2.35),	0, 	ZEROES,0,NULL, NULL};
planets sphinx2   ={ 3,	2.25,	0.02,	10.5,	EM(2.35),	FALSE,	EM(2.35),	0, 	ZEROES,0,NULL, &sphinx3};
planets sphinx    ={ 2,	1.6,	0.02,	10.5,	EM(2.2),	FALSE,	EM(2.2),	0, 	ZEROES,0,NULL, &sphinx2};
planets manticore ={ 1,	1.115,	0.017,	23.5,	EM(1.01),	FALSE,	EM(1.01),	0, 	ZEROES,0,NULL, &sphinx};


star	manticores[] = 
// L		Mass	Mass2	Eccen.	SMAxis	 Planets	Designation			Name
{{1.00,		1.00,	0,		0,		0,		 &mercury,	"Sol",		 	 1, "The Solar System"},
 {1.24,		1.047,	1.0,	0.05,	79.2,	 &manticore,"Manticore A",	 1, "Manticore A"},
 {1.0,		1.00,	1.047,	0.05,	79.2,	 NULL,		"Manticore B",	 1, "Manticore B"},
};

catalog	manticore_cat	= {sizeof(manticores) / sizeof (star),	"B", &manticores};

star	helios[] = 
// L		Mass	Mass2	Eccen.	SMAxis	 Planets	Designation		Name
{{1.00,		1.00,	0,		0,		0,		 &mercury,	"Sol",		 1, "The Solar System"},
 {1.08,		1.0,	0.87,	0.45,	8.85,	 NULL,		"Helio A",	 1, "Helio A"},
 {0.83,		0.87,	1.0,	0.45,	8.85,	 NULL,		"Helio B",	 1, "Helio B"},
};

catalog	helio		= {sizeof(helios) / sizeof (star), "?",	&helios};

			     /*	No Orbit Eccen. Tilt   Mass    Gas Giant? Dust Mass   Gas */
planets ilaqrb={1, 0.21, 0.1,   0,     EM(600.),TRUE,     0,   EM(600.), ZEROES,0,NULL, NULL};
planets ilaqrc={2, 0.13, 0.27,  0,     EM(178.),TRUE,     0,   EM(178.), ZEROES,0,NULL, &ilaqrb};
planets ilaqrd={3, 0.021,0.22,  0,     EM(5.9), FALSE,    EM(5.9),    0, ZEROES,0,NULL, &ilaqrc};	// EM(5.9) or 7.53 +/- 0.70 Earth-masses

star	ilAqrs[] = 
// L		Mass	Mass2	Eccen.	SMAxis	 Planets	Designation	Celes	Name
{{1.00,		1.00,	0,		0,		0,		 &mercury,	"Sol",		1, "The Solar System"},
{0.0016,	0.32,	0,		0,		0,		 &ilaqrd,	"IL Aqr",	1, "IL Aquarii/Gliese 876"}	// 15.2
};

catalog	ilAqr_cat		= {sizeof(ilAqrs) / sizeof (star),	"G", &ilAqrs};

//void usage(char*);

void usage(char *prognam)
{
// help - usage
	fprintf(stderr, "Usage: %s [options] [system name]\n", prognam);
	fprintf(stderr, "  Options:\n"
					"    -s#  Set random number seed\n"
					"    -m#  Specify stellar mass\n"
					"    -n#  Specify number of systems\n"
					"    -i#  Number to increment each seed by\n"
					"    -x   Use the Solar System's masses/orbits\n"
					"    -d   Use Dole's %d nearby stars\n"
					"    -D#  Use Dole's system #n\n"
					"    -w   Use %d nearby stars taken from the Web\n"
					"    -W#  Use Web system #n\n"
					"    -p   Specify the path for the output file\n"
					"    -o   Output file name\n"
					"    -t   Text-only output\n"
					"    -v#  Set verbosity (hex value)\n"
					"    -l   List nearby stars and exit\n"
					"    -H   Output only systems with habitable planets\n"
					"    -2   Only systems with 2 or more habitable planets\n"
					"    -E   Only systems with earthlike planets\n"
					"\n"
					"  Experimental options (may go away):\n"
					"    -c   Output Celestia .ssc file on stdout\n"
					"    -e   Output Excel .csv file\n"
					"    -V   Create vector graphics (SVG) system image\n"
					"    -k   Incorporate known planets (incomplete)\n"
					"          (use only orbital data at present)\n"
					"          Without -k, -c skips systems with known planets.\n"
					"    -g   Show atmospheric gases\n"
					"    -Z   Dump tables used for gases and exit\n"
					"    -M   Do moons (highly experimental)\n"
					"\n"
					"        Nearest stars taken from:\n"
					"          http://www.solstation.com/stars.htm\n"
                    "\n"
					"\n\n Additional long options in this version, under development:\n\n"
                    "    --help \n"
                    "    --mass <coeff> default: 1.0 mass of star\n"
                    "    --luminosity <coeff> default: 1.0 luminosity of star \n"
                    "    --gasdust <coeff> tex. --gasdust 50.0 gas per dust ratio \n Note: in long options you must precede space before value of parameter.\n"
                    "    --density <coeff> default: 1.0 density of dust, relative to default value\n"
                    "    --alpha <coeff> default: 5.0 dust density coeff alpha by Dole \n"
                    "    --gamma <coeff> default: 3.0 dust density distribution coeff gamma by Dole\n"
                    "    --beta <coeff> default: 1.2e-5 critical gas accretion mass  coefficient by Dole\n"
                    "    --disctype <int> default: 0, but can be 1 , 2 . \n With 1 you can set --alpha, with 2 you can set --alpha and --gamma , that are stdev and mean value of ring of surface density \n"
                    "    --nearest <coeff> default: 0.3 distence of nearest planet AU, before possibly migration\n"
                    "    --farthest <coeff> default: 50.0 farthest planet AU, unmigrated\n"
                    "    --discradius <coeff> default: 200.0 radius of dust disk AU\n"
                    "    --discecc <coeff> default: 0.2 eccentricity of cloud, greater : fwer planets\n"
                    "  \n"
                    "    --migrate <coeff> default: 1.0, that is no post-formation migration. \n tex --migrate 0.1 moves planets inward by coefficient 0.1"
                    "  \n"
                    "  Planet filtering out due to orbital disturbances by mutual gravitation:\n"
                    "  \n"
                    "   --filterhill <coeff>, default off \n"
                     "  --filterperiod <coeff> 2, default off "   
                    "\n"

					"\n"
					"        StarGen: %s\n"
					"\n",
					dole.count - 1,
					solstation.count - 1,
					stargen_revision);
}

int main (int argc, char *argv[])
{
	actions		action					= aGenerate;
	char		flag_char				= '?';
	char		path[300]				= SUBDIR;
	char 		url_path_arg[300]		= "";
	char		filename_arg[300]		= "";
	char 		arg_name [80] 			= "";

	int			use_stdout				= FALSE;
	char *		prognam;
	long double	mass_arg				= 0.0;
	long		seed_arg				= 0;
	int  		count_arg 				= 1;
	int			increment_arg			= 1;
	catalog	*	catalog					= NULL;
	int 		sys_no_arg				= 0;

	long double	ratio_arg				= 0.0;

	int			flags_arg				= 0;
	int			out_format				= ffHTML;
	int			graphic_format			= gfGIF;
	
	char 		*c						= NULL;
	int  		skip					= FALSE;
	int  		index					= 0;


    double dust_density_coeff2=1.0;
    double dabuffer=0.0;

#ifdef macintosh
	_ftype 		= 'TEXT';
	_fcreator 	= 'R*ch';
	argc = ccommand (&argv);
#endif
	
	prognam = argv[0];
	
	if ((c = strrchr(prognam, DIRSEP[0])) != NULL)
		prognam = c + 1;
	
	if (argc <= 1)
	{
		usage(prognam);
		return(1);
	}
	
	while (--argc > 0 && (*++argv)[0] == '-') {
		for (c = argv[0]+1, skip=FALSE; 
			 (*c != '\0') && (!(skip)); 
			 c++)
			switch (*c) 
			{
			//case '-':
			//	use_stdout = TRUE;
			//	break;
        case '-':  // "--" case


         if (strcmp(argv[0], "--gasdust") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                            //gasdust = atof(argv[1]);
                            sscanf (argv[1], "%lf", &dabuffer);
                            gasdust=(long double)dabuffer;
                            printf("gasdust %f", gasdust);
                            //exit(-1);

                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    }
            else if (strcmp(argv[0], "--density") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            dust_density_coeff2=(double)dabuffer;
                            printf(" dust_density_coeff2 %f",  dust_density_coeff2);
                            //double dust_density_coeff2=1.0;
                
                            dust_density_coeff=(long double)DUST_DENSITY_COEFF* (long double)dust_density_coeff2;
                            printf("\ndust_density_coeff %lf", dust_density_coeff);
                           // exit(-1);

                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 

            else if (strcmp(argv[0], "--migrate") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            migratek=(double)dabuffer;
                            printf(" migration coefficient %f",  migratek);
                            //double dust_density_coeff2=1.0;
                


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 

            else if (strcmp(argv[0], "--filterhill") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            BHILL_CRITERION =(double)dabuffer;
                            printf(" migration coefficient %f",  migratek);
                            //double dust_density_coeff2=1.0;
                
                            USE_FILTERING = 1;
                            FILTER_ASTEROIDS = 0; // Set to 1 to remove asteroids, 0 to keep them
                            USE_HILL = 1;         // Set to 1 to use Hill criterion, 0 for Orbital Perio

                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 
            else if (strcmp(argv[0], "--filterperiod") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            BHILL_CRITERION =(double)dabuffer;
                            printf(" migration coefficient %f",  migratek);
                            //double dust_density_coeff2=1.0;
                
                            USE_FILTERING = 1;
                            FILTER_ASTEROIDS = 0; // Set to 1 to remove asteroids, 0 to keep them
                            USE_HILL = 0;         // Set to 1 to use Hill criterion, 0 for Orbital Perio

                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 
  
            else if (strcmp(argv[0], "--alpha") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            alphaa=(long double)dabuffer;
                            printf(" alpha coeff  %f",  alphaa);
                            //double dust_density_coeff2=1.0;
                


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 

            else if (strcmp(argv[0], "--gamma") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            gammaa=(long double)dabuffer;
                            printf(" gamma coeff  %f", gammaa);
                            //double dust_density_coeff2=1.0;
                


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 
   
         else if (strcmp(argv[0], "--nearest") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            nearestk=(long double)dabuffer;
                            printf(" nearest planet  %f", nearestk);
                            //double dust_density_coeff2=1.0;
                


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 
         else if (strcmp(argv[0], "--farthest") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            farthestk=(long double)dabuffer;
                            printf(" farthest planet  %f", farthestk);
                            //double dust_density_coeff2=1.0;
                


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 
         else if (strcmp(argv[0], "--discradius") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            diskradiusk=(long double)dabuffer;
                            printf(" dust disk radius  %f", diskradiusk);
                            //double dust_density_coeff2=1.0;
                


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 
         else if (strcmp(argv[0], "--discecc") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            cloud_eccentricity=(long double)dabuffer;
                            printf(" disc cloud eccentricity %f", cloud_eccentricity);
                            //double dust_density_coeff2=1.0;
                


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 

         else if (strcmp(argv[0], "--beta") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            sscanf (argv[1], "%lf", &dabuffer);
                            bee=(long double)dabuffer;
                            printf(" beta - critical mass param %f", bee);
                            //double dust_density_coeff2=1.0;
                


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 

         else if (strcmp(argv[0], "--disctype") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            //sscanf (argv[1], "%li", &dabuffer);
                            //use_exponent_disk=(int)atoi(argv[1]);
                             sscanf (argv[1], "%li", &use_exponent_disk);

                            printf(" selected disc surface density profile  type %i", use_exponent_disk);
                            //double dust_density_coeff2=1.0;
                            //exit(-1);


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 
         else if (strcmp(argv[0], "--mass") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                    
                            //sscanf (argv[1], "%li", &dabuffer);
                            //use_exponent_disk=(int)atoi(argv[1]);
                             sscanf (argv[1], "%li", &mass_arg);

                            printf(" mass %i", mass_arg);
                            //double dust_density_coeff2=1.0;
                            //exit(-1);


                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }
                    } 
         else if (strcmp(argv[0], "--resonances") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                           
                            sscanf (argv[1], "%lf", &dabuffer);
                            //use_exponent_disk=(int)atoi(argv[1]);
                           //  sscanf (argv[1], "%li", &base_resonance);
                            base_resonance=(long double) dabuffer;

                            printf(" base resonance  %f", dabuffer);
                            //double dust_density_coeff2=1.0;
                            //exit(-1);

                                        order_to_resonances=1;
                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }

                        else
                        {
                           order_to_resonances=1;
                              skip=TRUE;
                        }
                    } 



         else if (strcmp(argv[0], "--luminosity") == 0) {
                        // Käsitellään '--gasdust 50' tyyppinen argumentti
                        if (argc > 1) {
                           
                            sscanf (argv[1], "%lf", &dabuffer);
                            //use_exponent_disk=(int)atoi(argv[1]);
                           //  sscanf (argv[1], "%li", &base_resonance);
                            par_luminosity=(long double) dabuffer;

                            printf(" luminosity  %f", dabuffer);
                            //double dust_density_coeff2=1.0;
                            //exit(-1);

                            use_own_luminosity=1;
                            argv++;  // Siirrytään seuraavaan argumenttiin
                            argc--;
                            skip=TRUE;
                        }

          
                    } 

         else if (strcmp(argv[0], "--help") == 0) {            
                            usage("stargen2");
                            skip=TRUE;
                       
                    } 
         else if (strcmp(argv[0], "--usage") == 0) {            
                            usage("stargen2");
                            skip=TRUE;
                       
                    } 
                   


           else {
                        use_stdout = 1;  // "--" ilman muuta tarkastelua
                    }      

			case 's':	// set random seed
				seed_arg = atol(&(*++c));
				skip = TRUE;
				break;
			case 'm':	// set mass of star
			{
				double	m;	// gnu C doesn't like to scanf long doubles
				
				sscanf (++c, "%lf", &m);
				mass_arg = m;
				
				skip = TRUE;
				break;
			}
		





		case 'n':	// number of systems
				count_arg = atoi(&(*++c));
				skip = TRUE;
				break;
			case 'i':	// number of systems
				increment_arg = atoi(&(*++c));
				skip = TRUE;
				break;
/*
			case 'T':	// Use the solar system with Titan, not Saturn
				jupiter.next_planet = &titan2;
 */
			case 'x':	// Use the solar system
				flag_char = *c;
				flags_arg |= fUseSolarsystem;
				if (mass_arg == 0.0)
					mass_arg = 1.0;
				break;
			case 'a':	// Use the solar system varying earth
				flag_char = *c;
				flags_arg |= fReuseSolarsystem;
				break;
			case 'D':
				catalog = &dole;
				flag_char = toupper(*c);
				++c;
				if ((toupper(*c) != 'X') && (*c != '\0'))
					sys_no_arg = atoi(c) + 1;
				skip = TRUE;
				break;
			case 'W':
				catalog = &solstation;
				flag_char = toupper(*c);
				++c;
				if ((toupper(*c) != 'X') && (*c != '\0'))
					sys_no_arg = atoi(c) + 1;
				skip = TRUE;
				break;
			case 'F':
				catalog = &jimb;
				flag_char = toupper(*c);
				++c;
				if ((toupper(*c) != 'X') && (*c != '\0'))
					sys_no_arg = atoi(c) + 1;
				skip = TRUE;
				break;
			case 'f':
				catalog = &jimb;
				flag_char = toupper(*c);
				break;
			case 'd':
				catalog = &dole;
				flag_char = toupper(*c);
				break;
			case 'w':
				catalog = &solstation;
				flag_char = toupper(*c);
				break;
			case 'b':						// experimental catal (Manticore, Helios etc.)
				catalog = &manticore_cat;
				flag_char = toupper(*c);
				break;
			case 'B':
				catalog = &manticore_cat;
				flag_char = toupper(*c);
				++c;
				if ((toupper(*c) != 'X') && (*c != '\0'))
					sys_no_arg = atoi(c) + 1;
				skip = TRUE;
				flags_arg |= fNoGenerate;
				sphinx.greenhouse_effect = TRUE;
				break;
			case 'G':
				catalog = &ilAqr_cat;
				flag_char = toupper(*c);
				++c;
				if ((toupper(*c) != 'X') && (*c != '\0'))
					sys_no_arg = atoi(c) + 1;
				skip = TRUE;
				flags_arg |= fNoGenerate;
				break;
			case 'o':
				if (*++c == '\0')
					if (argc > 1)
					{
						--argc;
						c = (++argv)[0];
					}
				
				if (*c != '\0')
					strcpy(filename_arg, c);

				skip = TRUE;
				break;
			case 't':	// display text
				out_format = ffTEXT;
				break;
			case 'e':
				out_format = ffCSV;
				break;
			case 'C':
				out_format = ffCSVdl;
				break;
			case 'c':
				out_format = ffCELESTIA;
				break;
			case 'V':
				graphic_format = gfSVG;
				break;
			case 'S':
				graphic_format = gfSVG;
				out_format = ffSVG;
				break;
			case 'k':
				flags_arg |= fUseKnownPlanets;
				break;
			case 'p':
				if (*++c == '\0')
					if (argc > 1)
					{
						--argc;
						c = (++argv)[0];
					}
				
				if (*c != '\0')
					strcpy(path, c);
				
				if (strcmp(path + strlen(path) - strlen(DIRSEP), DIRSEP) != 0)
					strncat (path, DIRSEP, 80-strlen(path));
					
				skip = TRUE;
				break;
			case 'u':
				if (*++c == '\0')
					if (argc > 1)
					{
						--argc;
						c = (++argv)[0];
					}
				
				if (*c != '\0')
					strcpy(url_path_arg, c);
				
				if (strcmp(url_path_arg + strlen(url_path_arg) - strlen("/"), "/") != 0)
					strncat (url_path_arg, "/", 80-strlen(url_path_arg));
				
				skip = TRUE;
				break;
			case 'g':
				flags_arg |= fDoGases;
				break;
			case 'v':	// verbosity
				if (isdigit(*(c+1)))
				{
					sscanf (++c, "%x", &flag_verbose);
					skip = TRUE;
					if (flag_verbose & 0x0001)
						flags_arg |= fDoGases;
				}
				else
					action = aListVerbosity;
				break;
			case 'l':
				action = aListCatalog;
				break;
			case 'L':
				action = aListCatalogAsHTML;
				break;
			case 'z':
				action = aSizeCheck;
				break;
			case 'Z':
				action = aListGases;
				break;
			case 'M':
				flags_arg |= fDoMoons;
				break;
			case 'H':
				flags_arg |= fDoGases | fOnlyHabitable;
				break;
			case '2':
				flags_arg |= fDoGases | fOnlyMultiHabitable;
				break;
			case 'J':
				flags_arg |= fDoGases | fOnlyJovianHabitable;
				break;
			case 'E':
				flags_arg |= fDoGases | fOnlyEarthlike;
				break;
			case 'A':
			{
				double ratio;
				
				sscanf (++c, "%lf", &ratio);
				skip = TRUE;
				
				if (ratio > 0.0)
					ratio_arg = ratio;
				break;
			}
			
			default:
				fprintf (stderr, "Unknown option: %s\n", c);
			case '?':
			case 'h':
				usage(prognam);
				return (1);
			}
	}
	
	for (index = 0; index < argc; index++) {
		if ((strlen(argv[index]) + strlen(arg_name)) < sizeof(arg_name))
		{
			if (strlen(arg_name))
				strcpy(arg_name+strlen(arg_name), " ");
			
			strcpy(arg_name+strlen(arg_name), argv[index]);
		}
	}
	
 // printf("%lf", (double)mass_arg);

 //   exit(-1);

	stargen (action,
			 flag_char,
			 path,
			 url_path_arg,
			 filename_arg,
			 arg_name,
			 
			 use_stdout ? stdout : NULL,
			 stderr,
			 prognam,
			 mass_arg,
			 seed_arg,
			 count_arg,
			 increment_arg,
			 catalog,
			 sys_no_arg,
			 
			 ratio_arg,
			 
			 flags_arg,
			 out_format,
			 graphic_format
			 );

	return(0);
}
