#ifndef FC_HEADER_INCLUDED
#define FC_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define FC_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define FC_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define FC_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/* Mangling for Fortran module symbols with underscores. */
#define FC_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define curfit FC_GLOBAL(curfit, CURFIT)
#define splev FC_GLOBAL(splev, SPLEV)
#define splder FC_GLOBAL(splder, SPLDER)
#define surfit FC_GLOBAL(surfit, SURFIT)
#define bispev FC_GLOBAL(bispev, BISPEV)
#define parder FC_GLOBAL(parder, PARDER)

#endif
