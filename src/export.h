/* export.h */

/*
 *  Copyright (c) 1992-1998, Research Systems Inc.  All rights reserved.
 *  Reproduction by any means whatsoever is prohibited without express
 *  written permission.
*/

/*
 * NOTICE OF INTENT. PLEASE READ:
 *
 *     This file is intended to supply the external definitions
 *     required to write and link code against IDL, either to extend
 *     IDL (CALL_EXTERNAL, Dynamically loadable modules, LINKIMAGE) or
 *     to use IDL to extend your own program (Callable IDL). Our goal is
 *     that including this file in your program will enable your program
 *     to access the public interfaces described in the IDL documentation,
 *     especially the External Development Guide.
 *
 *     NOT EVERYTHING IN THIS FILE IS A STABLE PUBLIC INTERFACE.
 *     For various technical reasons, there are things in this file that
 *     are considered to be private to Research Systems and subject to
 *     immediate change without notice. Anything in this file that is not
 *     explicitly documented elsewhere falls into this category. You should
 *     not use such interfaces, as they are not publically supported and
 *     may change or disappear without warning.
 *
 *     Examples of such interfaces are:
 *         - Types and constants that are not discussed in the documentation.
 *         - Functions that are not discussed elsewhere.
 *         - Fields of structures that are not discussed in the documentation
 *           even though the structure itself *is* documented.
 */

#ifndef export_IDL_DEF
#define export_IDL_DEF

#ifdef __cplusplus
    extern "C" {
#endif


/***** Definitions from msg_code *****/

#ifndef msg_code_IDL_DEF
#define msg_code_IDL_DEF

/* Warning: These codes can change between releases. */

#define IDL_M_GENERIC                    -1
#define IDL_M_NAMED_GENERIC              -2
#define IDL_M_SYSERR                     -4
#define IDL_M_NOSUPPORT                  -5
#define IDL_M_NOPROCSUPPORT              -6

#endif                         /* msg_code_IDL_DEF */




/***** Definitions from config *****/

#ifndef config_IDL_DEF
#define config_IDL_DEF

#ifndef IDL_DEBUGGING
#define IDL_DEBUGGING 1
#endif


#ifdef sun
#ifdef SUN_OS_4_1
#define SUN_BSD
#else				/* SunOS 5.0 and later is SVR4 based */
#define SUN_SYSV
#ifdef sparc
#define SUN_SYSV_SPARC
#else
#define SUN_SYSV_X86
#endif
#endif
#endif				/* sun */

#ifdef vms
#ifdef vax
#define VAX_VMS
#else				/* vax */
#define ALPHA_VMS
#ifdef __IEEE_FLOAT
#define VMS_IEEE		/* ALPHA/VMS with IEEE arithmetic */
#endif
#endif				/* vax */
#endif				/* vms */

#if defined(__alpha) && defined(__osf__)
#define ALPHA_OSF
#endif

#if defined(_ALPHA_) && defined(WIN32)
#define ALPHA_WIN32
#endif

/*
 * Proper ANSI C doesn't allow pre-defined cpp macros that don't start
 * with an underscore. Explictly define the ones we use.
 */
#if !defined(unix) && (defined(__unix__) || defined(__unix) || defined(_AIX))
#define unix
#endif
#if defined(__hpux) && !defined(hpux)
#define hpux
#endif
#if defined(__hp9000s800) && !defined(hp9000s800)
#define hp9000s800
#endif

/* MAC is not predefined by Macintosh compilers but is used widely by IDL */
#if defined(__MACOS__) && !defined(MAC)
#define MAC
#endif


/*
 * In IDL terms, a longword is 32 bits. By using this typedef in place
 * of the standard C "long", we can handle systems where this assumption
 * does not hold.
 */
#ifdef ALPHA_OSF
#define LONG_NOT_32		/* Can be used to configure code */
#endif


/*
 * IDL is built with ANSI C on all platforms. However, there are
 * still a few K&R compilers out there. If you are using such a
 * compiler, define the preprocessor symbol IDL_CC_NOT_ANSI before
 * including this header file.
 *
 * The following macro will suppress argument prototypes in function
 * declarations allowing them to compile with a K&R compiler. Please be
 * aware, however, that the actual functions were compiled with an ANSI
 * compiler with full prototypes, so the default K&R type promotion rules
 * will cause runtime errors if you pass char, short, or float arguments.
 * K&R C widens these to int or double, but the actual routines expect
 * the un-widened types. Routines with these argument types should not
 * be called from K&R compiled code.
 */
#ifdef IDL_CC_NOT_ANSI
#define IDL_ARG_PROTO(args) ()
#else
#define IDL_ARG_PROTO(args) args
#endif


#ifdef FALSE
#undef FALSE
#endif
#define FALSE   (0)

#ifdef TRUE
#undef TRUE
#endif
#define TRUE    (1)

/*
 * The following definitions are to be used in all modules. They
 * are used by cx to generate .x files. cx_public is used for functions that
 * are used outside a module, but which are RSI private. cx_export is used
 * for functions that people outside of RSI can use.
*/

#define cx_public		/* C default file scope, private to RSI */
#define cx_export		/* C default file scope, visible in export.h */

/*
 * Microsoft Windows has two calling conventions, __stdcall and __cdecl. We
 * refer to them using the following macros so that they can be easily
 * nulled on non-MS platforms.
 */
#ifdef WIN32
#define IDL_STDCALL __stdcall
#define IDL_CDECL __cdecl
#else
#define IDL_STDCALL
#define IDL_CDECL
#endif
 

/*
 * Some platforms support tty based user interaction (Unix, VMS)
 * while others (Macintosh, Windows) don't. On those that don't, we
 * want to avoid compiling code that will never be called.
 * This symbol is defined on those systems that support ttys.
 */
#if defined(unix) || defined(VMS)
#define IDL_OS_HAS_TTYS
#endif



#if IDL_DEBUGGING < 2
#define IDL_REGISTER register	/* Use explicit register declarations */
#else
#define IDL_REGISTER
#endif


/**** Maximum # of array dimensions ****/
#define IDL_MAX_ARRAY_DIM 8

/**** Maximum # of params allowed in a call ****
  NEVER make this > 127.  */
#define IDL_MAXPARAMS 64

/**** Longest allowed file path specification ****/
#if defined(MAC) || defined(WIN32)
#define IDL_MAXPATH 255
#endif
#ifdef unix
#define IDL_MAXPATH    1024		/* That's what BSD allows */
#endif
#ifdef VMS
#define IDL_MAXPATH    264
#endif



#endif				/* config_IDL_DEF */




/***** Definitions from defs *****/

#ifndef defs_IDL_DEF
#define defs_IDL_DEF

#if !defined(WIN32) || !defined(PLTYPES)
typedef unsigned char UCHAR;	/* Unsigned character type */
#endif

/* Boolean. */ 
typedef enum {
    IDL_FALSE = 0,
    IDL_TRUE = 1
} IDLBool_t;

/*
 * IDL integer types. For historical reasons, the following types
 * do not have a typedef in our code and instead use the C primatives
 * listed here:
 *
 *	TYP_BYTE	UCHAR
 *	TYP_INT		short
 *
 * TYP_LONG used to be represented by "long", but the advent of 64-bit
 * systems where long can be 64-bits caused us to define a typedef for
 * it. All the newer types have typedefs as well.
 */
typedef unsigned short IDL_UINT;
#ifdef LONG_NOT_32
typedef int IDL_LONG;
typedef unsigned int IDL_ULONG;
#else
typedef long IDL_LONG;
typedef unsigned long IDL_ULONG;
#endif

#ifdef WIN32
typedef __int64 IDL_LONG64;
typedef unsigned __int64 IDL_ULONG64;
#else
typedef long long IDL_LONG64;
typedef unsigned long long IDL_ULONG64;
#endif

/* Type used for pointer and object reference variables is same as IDL_LONG */
typedef IDL_ULONG IDL_HVID;


/*
 * Define IDL_VARIABLE type values - Note that IDL_TYP_UNDEF is always 0 by
 * definition. It is correct to use the value 0 in place of IDL_TYP_UNDEF.
 * It is not correct to assume the value assigned to any other
 * type - the preprocessor definitions below must be used.
 */

#define IDL_TYP_UNDEF		0
#define IDL_TYP_BYTE            1
#define IDL_TYP_INT             2
#define IDL_TYP_LONG            3
#define IDL_TYP_FLOAT           4
#define IDL_TYP_DOUBLE          5
#define IDL_TYP_COMPLEX         6
#define IDL_TYP_STRING          7
#define IDL_TYP_STRUCT          8
#define IDL_TYP_DCOMPLEX        9
#define IDL_TYP_PTR		10
#define IDL_TYP_OBJREF		11
#define IDL_TYP_UINT		12
#define IDL_TYP_ULONG		13
#define IDL_TYP_LONG64		14
#define IDL_TYP_ULONG64		15


#define IDL_MAX_TYPE            15
#define IDL_NUM_TYPES           16

/*
 * Various machines use different data types for representing memory
 * and file offsets and sizes. We map these to IDL types using the
 * following definitions. Doing it this way lets us easily change
 * the mapping here without having to touch all the code that uses
 * these types.
 */

/*
 * Memory is currently limited to 2^31 on all platforms. If this should
 * be increased to 64-bits, IDL_MEMINT_64 would be defined to allow code
 * to work around it.
 *
 * MEMINT must always be a signed type.
 */
#define IDL_TYP_MEMINT	 IDL_TYP_LONG
#define IDL_TYP_UMEMINT	 IDL_TYP_ULONG
#define IDL_MEMINT	 IDL_LONG
#define IDL_UMEMINT	 IDL_ULONG

#if defined(SUN_SYSV) || defined(ALPHA_OSF) || defined(sgi) || defined(hpux)
				/* Files can have 64-bit sizes */
#define IDL_FILEINT_64
#define IDL_TYP_FILEINT	  IDL_TYP_LONG64
#define IDL_FILEINT	  IDL_LONG64
#else				/* Stick with 2^31 sized files */
#define IDL_TYP_FILEINT	  IDL_TYP_LONG
#define IDL_FILEINT	  IDL_LONG
#endif




/*
 * The above type codes each have a bit mask value associated with
 * them. The bit mask value is computed as (2**Type_code), but the
 * following definitions can also be used. Some routines request the
 * bit mask value instead of the type code value.
 *
 * Simple types are everything except TYP_STRUCT, TYP_PTR, and TYP_OBJREF.
*/
#define IDL_TYP_B_SIMPLE            62207
#define IDL_TYP_B_ALL               65535

/* This macro turns it's argument into its bit mask equivalent.
 * The argument type_code should be one of the type codes defined
 * above.
*/ 
#define IDL_TYP_MASK(dim_code)      (1 << dim_code)



/***** IDL_VARIABLE flag values ********/

#define IDL_V_CONST         1
#define IDL_V_TEMP          2
#define IDL_V_ARR           4
#define IDL_V_FILE          8
#define IDL_V_DYNAMIC       16
#define IDL_V_STRUCT        32
#define IDL_V_NOT_SCALAR    (IDL_V_ARR | IDL_V_FILE | IDL_V_STRUCT)

/**** IDL_ARRAY flag values ****/
#define IDL_A_FILE          1	/* Array is a FILE variable (ASSOC) */
#define IDL_A_NO_GUARD      2	/* Indicates no data guards for array */
#define IDL_A_FILE_PACKED   4	/* If array is a FILE variable and the data
				   type is IDL_TYP_STRUCT, then I/O to
				   this struct should assume packed data
				   layout compatible with WRITEU instead of
				   being a direct mapping onto the struct
				   including its alignment holes. */


/**** Basic IDL structures: ****/

typedef struct {
  float r,i;
} IDL_COMPLEX;

typedef struct {
  double r,i;
} IDL_DCOMPLEX;

typedef struct {		/* Define string descriptor */
  unsigned short slen;		/* Length of string, 0 for null */
  short stype;			/* type of string, static or dynamic */
  char *s;			/* Addr of string */
} IDL_STRING;


/**** IDL identifiers ****/
typedef struct _idl_ident {
  struct _idl_ident *hash;	/* Must be the first field */
  char *name;                   /* Identifier text (NULL terminated */
  int len;			/* # of characters in id, not counting NULL
				   termination. */
} IDL_IDENT;


/*
 * Type of the free_cb field of IDL_ARRAY. When IDL deletes a variable and
 * the free_cb field of ARRAY non-NULL, IDL calls the function that field
 * references, passing the value of the data field as it's sole argument.
 * The primary use for this notification is to let programs know when
 * to clean up after calls to IDL_ImportArray(), which is used to create
 * arrays using memory that IDL does not allocate.
 */
typedef void (* IDL_ARRAY_FREE_CB) IDL_ARG_PROTO((UCHAR *));

/* Type of the dim field of an IDL_ARRAY. */
typedef IDL_MEMINT IDL_ARRAY_DIM[IDL_MAX_ARRAY_DIM];

typedef struct {		/* Its important that this block
				   be an integer number of longwords
				   in length to ensure that array
				   data is longword aligned.  */
  IDL_MEMINT elt_len;		/* Length of element in char units */
  IDL_MEMINT arr_len;		/* Length of entire array (char) */
  IDL_MEMINT n_elts;		/* total # of elements */
  UCHAR *data;			/* ^ to beginning of array data */
  UCHAR n_dim;			/* # of dimensions used by array */
  UCHAR flags;			/* Array block flags */
  short file_unit;		/* # of assoc file if file var */
  IDL_ARRAY_DIM dim;		/* dimensions */
  IDL_ARRAY_FREE_CB free_cb;	/* Free callback */
  IDL_FILEINT offset;		/* Offset to base of data for file var */
  IDL_LONG data_guard;		/* Guard longword */
} IDL_ARRAY;

typedef struct {		/* Reference to a structure */
  IDL_ARRAY *arr;		/* ^ to array block containing data */
  struct _idl_structure *sdef;	/* ^ to structure definition */
} IDL_SREF;

/* IDL_ALLTYPES can be used to represent all IDL_VARIABLE types */
typedef union {
  char sc;			/* A standard char, where "standard" is defined
				   by the compiler. This isn't an IDL data
				   type, but having this field is sometimes
				   useful for internal code */
  UCHAR c;			/* Byte value */
  short i;			/* Integer short value */
  IDL_UINT ui;			/* Unsigned integer short value */
  IDL_LONG l;			/* Long value */
  IDL_ULONG ul;			/* Unsigned long value */
  IDL_LONG64 l64;		/* 64-bit integer value */
  IDL_ULONG64 ul64;		/* Unsigned 64-bit integer value */
  float f;			/* Floating value */
  double d;			/* Double value */
  IDL_COMPLEX cmp;		/* Complex value */
  IDL_DCOMPLEX dcmp;		/* Double complex value */
  IDL_STRING str;		/* String descriptor */
  IDL_ARRAY *arr;		/* ^ to array descriptor */
  IDL_SREF s;			/* Descriptor of structure */
  IDL_HVID hvid;		/* Heap variable identifier */
  IDL_MEMINT memint;		/* Memory size or offset */
  IDL_FILEINT fileint;		/* File size or offset */
} IDL_ALLTYPES;

typedef struct {		/* IDL_VARIABLE definition */
  UCHAR type;			/* Type byte */
  UCHAR flags;			/* Flags byte */
  IDL_ALLTYPES value;
} IDL_VARIABLE;
typedef IDL_VARIABLE *IDL_VPTR;

typedef void (* IDL_PRO_PTR)();	/* ^ to interpreter procedure (ret is void) */
				/* ^ to interp. function (ret ^ to VAR) */
typedef IDL_VARIABLE *(* IDL_FUN_RET)();

/* Possible values for the flags field of IDL_SYSFUN_DEF */
#define IDL_SYSFUN_DEF_F_OBSOLETE	1   /* Routine is obsolete */
#define IDL_SYSFUN_DEF_F_KEYWORDS	2   /* Routine accepts keywords */
#define IDL_SYSFUN_DEF_F_METHOD		32   /* Routine is an object method */

/* This structure defines the format of a system procedure
   or function table entry. */
typedef struct {		/* System function definition */
  IDL_FUN_RET funct_addr;	/* Address of function, or procedure.  */
  char *name;			/* The name of the function */
  UCHAR arg_min;		/* Minimum allowed argument count. */
  UCHAR arg_max;		/* Maximum argument count.  The top
				   bit in arg_min is set to indicate that
				   the routine accepts keywords. */
  UCHAR flags;			/* IDL_SYSFUN_DEF_F_* flags */
} IDL_SYSFUN_DEF;

/*
 * Setting the top bit in the arg_min field of an IDL_SYSFUN_DEF passed to
 * IDL_AddSystemRoutine() is equivalent to the setting the
 * IDL_SYSFUNDEF_F_KEYWORDS bit in the flags field. This is strictly for
 * backwards compatability. Direct use of the flags field is preferred.
 */
#define IDL_KW_ARGS 128		/* Bit set in argmin indicating kw's allowed */


/*
 * Type of pointer to an IDL structure definition. This is an opaque type
 * not to be directly referenced.
 */
typedef void *IDL_StructDefPtr;

#endif				/* defs_IDL_DEF */




/***** Definitions from message *****/

#ifndef message_IDL_DEF
#define message_IDL_DEF

/*
 * action parameter to message is composed of two parts. These two masks
 * are used to separate them out.
 */
#define IDL_MSG_ACTION_CODE	0x0000ffff
#define IDL_MSG_ACTION_ATTR	0xffff0000

/* Allowed codes for action parameter to IDL_Message() */
#define IDL_MSG_RET	    0   /* Return to caller */
#define IDL_MSG_EXIT	    1   /* Terminate process via exit(3) */
#define IDL_MSG_LONGJMP	    2   /* General error. Obey the error handling
				   established by the ON_ERROR user
				   procedure. */
#define IDL_MSG_IO_LONGJMP  3   /* I/O error. Obey the error handling
				   established by the ON_IOERROR user
				   procedure. */
#define IDL_MSG_INFO	    4   /* Informational. Like IDL_MSG_RET, but won't
				   set !ERR or !ERR_STRING. Also,
				   inhibited by !QUIET */


/* Allowed attribute masks that can be OR'd into the action code */
#define IDL_MSG_ATTR_NOPRINT  0x00010000   /* Suppress the printing of
					      the error text to stderr,
					      but do everything else in
					      the normal way. */
#define IDL_MSG_ATTR_MORE  0x00020000   /* Use IDL_more() instead of printf(3S)
					   to output the message. The calling
					   routine must worry about calling
					   IDL_more_reset(). A side effect of
					   this is that the message goes to the
					   file named in IDL_more_reset(), not
					   necessarily stderr. */
#define IDL_MSG_ATTR_NOPREFIX  0x00040000   /* Don't output the normal
					       message (from
					       !ERROR_STATE.MSG_PREFIX),
					       just the message text. */
#define IDL_MSG_ATTR_QUIET     0x00080000   /* If the message would normally
					       be printed and !QUIET is
					       non-zero, the printing
					       is suppressed. Everything
					       else is updated as expected. */
#define IDL_MSG_ATTR_NOTRACE   0x00100000   /* Suppress traceback message */
#define IDL_MSG_ATTR_BELL      0x00200000   /* Ring the bell */
#define IDL_MSG_ATTR_SYS       0x00400000   /* IDL_Message() only: Include
					       system supplied error message,
					       via errno or whatever source
					       is appropriate for the
					       current platform. */

/* The type of elements in the defs argument to IDL_MessageDefineBlock() */
typedef struct {
  char *name;
  char *format;
} IDL_MSG_DEF;

/* Type returned by IDL_MessageDefineBlock() */
typedef void *IDL_MSG_BLOCK;


#endif				/* message_IDL_DEF */




/***** Definitions from macros *****/

#ifndef macros_IDL_DEF
#define macros_IDL_DEF

/* General math macros */
#define IDL_MIN(x,y) (((x) < (y)) ? (x) : (y))
#define IDL_MAX(x,y) (((x) > (y)) ? (x) : (y))
#define IDL_ABS(x) (((x) >= 0) ? (x) : -(x))

/* Return x in the range of min <= x <= max */
#define IDL_CLIP_TO_RANGE(x, min, max) \
  ((x) < (min) ? (min) : ((x) > (max) ? (max) : (x)))

/* Round x up modulo m. m must be a power of 2 : */
#define IDL_ROUND_UP(x,m) \
  (((x) + (m-1)) & (~(m-1)))

/**** Cast a pointer to (char *) ****/
#define IDL_CHAR(x) ((char *) x)

/**** Take the address of a variable and cast to a desired type ****/
#define IDL_CHARA(x) ((char *) &(x))
#define IDL_UCHARA(x) ((UCHAR *) &(x))
#define IDL_SHORTA(x) ((short *) &(x))
#define IDL_INTA(x) ((int *) &(x))
#define IDL_LONGA(x) ((IDL_LONG *) &(x))

/**** Get pointer to a valid string from an IDL_STRING descriptor */
#define IDL_STRING_STR(desc) ((desc)->slen ? (desc)->s : "")

#define IDL_DELTMP(v) { if ((v->flags) & IDL_V_TEMP) IDL_Deltmp(v); }

/**** How many elements are there in a C 1D array? ****/
#define IDL_CARRAY_ELTS(arr) (sizeof(arr)/sizeof(arr[0]))

#define IDL_EXCLUDE_UNDEF(v) { if (!v->type) \
	IDL_MessageVE_UNDEFVAR(v, IDL_MSG_LONGJMP); }
#define IDL_EXCLUDE_CONST(v)        { if (v->flags & IDL_V_CONST) \
	IDL_MessageVE_NOCONST(v, IDL_MSG_LONGJMP); }
#define IDL_EXCLUDE_EXPR(v)  { if (v->flags & (IDL_V_CONST | IDL_V_TEMP)) \
	IDL_MessageVE_NOEXPR(v, IDL_MSG_LONGJMP); }
#define IDL_EXCLUDE_FILE(v) { if (v->flags & IDL_V_FILE) \
	IDL_MessageVE_NOFILE(v, IDL_MSG_LONGJMP); }
#define IDL_EXCLUDE_STRUCT(v)       { if (v->flags & IDL_V_STRUCT) \
	IDL_MessageVE_NOSTRUCT(v, IDL_MSG_LONGJMP); }
#define IDL_EXCLUDE_FILE_OR_STRUCT(v) {if(v->flags &(IDL_V_FILE|IDL_V_STRUCT))\
	IDL_VarExclude(v, IDL_TYP_MASK(TYP_STRUCT), FALSE, FALSE, TRUE);}
#define IDL_EXCLUDE_COMPLEX(v)      { if ((v->type == IDL_TYP_COMPLEX) \
					  || (v->type == IDL_TYP_DCOMPLEX)) \
	IDL_MessageVE_NOCOMPLEX(v, IDL_MSG_LONGJMP); }
#define IDL_EXCLUDE_STRING(v)       { if (v->type == IDL_TYP_STRING)  \
	IDL_MessageVE_NOSTRING(v, IDL_MSG_LONGJMP); }
#define IDL_EXCLUDE_SCALAR(v) { if (!(v->flags & IDL_V_NOT_SCALAR)) \
	IDL_MessageVE_NOSCALAR(v, IDL_MSG_LONGJMP);}


/**** Ensure that variables possess certain attributes ****/
#define IDL_ENSURE_ARRAY(v) { if (!(v->flags & IDL_V_ARR)) \
	IDL_MessageVE_NOTARRAY(v, IDL_MSG_LONGJMP); }
#define IDL_ENSURE_SCALAR(v) { if (v->flags & IDL_V_NOT_SCALAR) \
	IDL_MessageVE_NOTSCALAR(v, IDL_MSG_LONGJMP);}
#define IDL_ENSURE_STRING(v) { if (v->type != IDL_TYP_STRING)  \
	IDL_MessageVE_REQSTR(v, IDL_MSG_LONGJMP);}
#ifndef IDL_ENSURE_SIMPLE
#define IDL_ENSURE_SIMPLE(v) IDL_VarEnsureSimple(v)
#endif
#define IDL_ENSURE_STRUCTURE(v) { if (!(v->flags & IDL_V_STRUCT)) \
	IDL_MessageVE_STRUC_REQ(v, IDL_MSG_LONGJMP);}
#define IDL_ENSURE_PTR(v) { if (v->type != IDL_TYP_PTR)  \
	IDL_MessageVE_REQPTR(v, IDL_MSG_LONGJMP);}
#define IDL_ENSURE_OBJREF(v) { if (v->type != IDL_TYP_OBJREF)  \
	IDL_MessageVE_REQOBJREF(v, IDL_MSG_LONGJMP);}


     /* Check if var has a dynamic part. If so, delete it using IDL_Delvar  */
#define IDL_DELVAR(v) { if (((v)->flags) & IDL_V_DYNAMIC) IDL_Delvar(v); }


/*
 * For systems that lack the bcopy() routines, mimic them using
 * the memcpy() routines. Under VMS, DEC C will inline these
 * because string.h is included.
 */
#if defined(SUN_SYSV) || !defined(unix)
#include <string.h>
#ifndef bcopy
#define bcopy(src,dest,len)     (memcpy((dest), (src), (len)))
#endif
#ifndef bzero
#define bzero(dest,len)         (memset((dest), 0, (len)))
#endif
#ifndef bcmp
#define bcmp(b1,b2,len)         (memcmp((b1), (b2), (len)))
#endif
#endif

#endif				/* macros_IDL_DEF */




/***** Definitions from crearr *****/

#ifndef crearr_IDL_DEF
#define crearr_IDL_DEF

/* The following define the valid values for the init arg to basic_array */
#define IDL_ARR_INI_ZERO   0	/* Zero data area */
#define IDL_ARR_INI_NOP    1	/* Don't do anything to data area */
#define IDL_ARR_INI_INDEX  2	/* Put 1-D index into each elt. */
#define IDL_ARR_INI_TEST   3	/* Test if enough memory is available */

/* Old names for the array init arg */
#define IDL_BARR_INI_ZERO  IDL_ARR_INI_ZERO
#define IDL_BARR_INI_NOP   IDL_ARR_INI_NOP
#define IDL_BARR_INI_INDEX IDL_ARR_INI_INDEX
#define IDL_BARR_INI_TEST  IDL_ARR_INI_TEST

#endif				/* crearr_IDL_DEF */




/***** Definitions from ez *****/

#ifndef ez_IDL_DEF
#define ez_IDL_DEF

/* These constants can be ORd together to form the value for the
   access field of IDL_EZ_ARG */
#define IDL_EZ_ACCESS_R     1	/* Arg is readable */
#define IDL_EZ_ACCESS_W     2	/* Arg is writable */
#define IDL_EZ_ACCESS_RW    3	/* Arg is readable and writable */


/* This macro turns it's argument into a bit mask suitable for
 * the allowed_dims field of IDL_EZ_ARG. The argument dim_code should be
 * 0 for scalar, 1 for 1D, 2 for 2D, etc...
 */ 
#define IDL_EZ_DIM_MASK(dim_code)   (1 << dim_code)

/* Define type mask of all numeric types: */
#define IDL_EZ_TYP_NUMERIC \
	( IDL_TYP_MASK(TYP_INT) | IDL_TYP_MASK(TYP_LONG) \
	 | IDL_TYP_MASK(TYP_FLOAT) | IDL_TYP_MASK(TYP_DOUBLE) \
	 | IDL_TYP_MASK(TYP_COMPLEX) | IDL_TYP_MASK(TYP_BYTE) \
	 | IDL_TYP_MASK(TYP_DCOMPLEX) | IDL_TYP_MASK(TYP_UINT) \
	 | IDL_TYP_MASK(TYP_ULONG) | IDL_TYP_MASK(TYP_LONG64) \
	 | IDL_TYP_MASK(TYP_ULONG64) )


/* These constants should be used instead of IDL_EZ_DIM_MASK when appropriate */
#define IDL_EZ_DIM_ARRAY    510	  /* Allow all but scalar */
#define IDL_EZ_DIM_ANY      511	  /* Allow anything */


/* These constants can be ORd together to form the value for the
   pre field of IDL_EZ_ARG. These actions are taken only if the argument
   has IDL_EZ_ACCESS_R. */
#define IDL_EZ_PRE_SQMATRIX         1	/* Arg must be a square matrix. */
#define IDL_EZ_PRE_TRANSPOSE        2	/* Transpose arg. This only happens
					   with read access. */


/* These constants can be ORd together to form the value for the
   post field of IDL_EZ_ARG. These actions are taken only if the argument
   has IDL_EZ_ACCESS_W.  If IDL_EZ_POST_WRITEBACK is not present, none of the
   other actions are considered, since that would imply wasted effort. */
#define IDL_EZ_POST_WRITEBACK       1	/* Transfer the contents of
					   uargv to the actual argument. */
#define IDL_EZ_POST_TRANSPOSE       2	/* Transpose uargv prior to writing. */


/*
 * IDL_EZ_ARG is the definition for the structure used by IDL_EzCall()
 * and IDL_EzCallCleanup() to define the plain arguments being passed
 * to a routine.
 */
typedef struct {
  short allowed_dims;		/* A bit mask that specifies the
				   allowed dimensions. Bit 0 means scalar,
				   bit 1 is 1D, etc. Use the EZ_DIM_* constants
				   defined in this file to specify this
				   value. */
  int allowed_types;		/* This is a bit mask defining the
				   allowed data types for the argument.
				   To convert the TYP_* type codes defined
				   in defs.h to the appropriate bits,
				   use the formula 2**(type_code) or use
				   the TYP_B_* bit masks defined in defs.h
				   NOTE: If you specify a value for convert,
				   its a good idea to specify IDL_TYP_B_ALL or
				   IDL_TYP_B_SIMPLE here. The type conversion
				   will catch any problems and your routine
				   will be more flexible. */
  short access;			/* Some combination of the EZ_ACCESS
				   constants defined above. */
  short convert;		/* If non-zero, the TYP_* type code to
				   which the argument will be converted.
				   A value of zero means that no conversion
				   will be applied. */
  short pre;			/* A bit mask that specifies special purpose
				   processing that should be performed on
				   the variable by IDL_EzCall(). These bits
				   are specified with the IDL_EZ_PRE_*
				   constants. This processing occurs *AFTER*
				   any type conversions specified by
				   convert. */
  short post;			/* A bit mask that specifies special purpose
				   processing that should be performed on
				   the variable by IDL_EzCallCleanup(). These
				   bits are specified with the IDL_EZ_POST_*
				   constants. */
  IDL_VPTR to_delete;		/* RESERVED TO EZ MODULE. DO NOT MAKE
				   USE OF OR CHANGE THIS FIELD. If EZ
				   allocated a temporary variable to satisfy
				   the conversion requirements given by the
				   convert field, the IDL_VPTR to that temp
				   is stashed here by IDL_EzCall for use by
				   IDL_EzCallCleanup(). */
  IDL_VPTR uargv;		/* After calling IDL_EzCall(), uargv contains
				   a pointer to the IDL_VARIABLE which is
				   the argument. */
  IDL_ALLTYPES value;		/* This is a copy of the value field
				   of the variable pointed at by uargv.
				   For scalar variables, it contains the
				   value, for arrays it points at the
				   array block. */
} IDL_EZ_ARG;

#endif				/* ez_IDL_DEF */




/***** Definitions from graphics *****/

#ifndef graphics_IDL_DEF
#define graphics_IDL_DEF

/* *** Structure defining current device and parameters: *** */
#define IDL_MAX_TICKN 60	/* Max # of axis annotations */

#define IDL_COLOR_MAP_SIZE 256	 /* Size of internal color map. */

#define IDL_NUM_LINESTYLES 6	/* # of line styles */
#define IDL_X0 0		/* Subscripts of fields for rect structures */
#define IDL_Y0 1
#define IDL_X1 2
#define IDL_Y1 3
#define IDL_Z0 4
#define IDL_Z1 5


#define IDL_AX_LOG 1		/* Axis type values */
#define IDL_AX_MAP 2		/* Old style maps */
#define IDL_AX_MAP1 3		/* New style maps */

#define IDL_AX_EXACT 1		/* Axis style values: */
#define IDL_AX_EXTEND 2
#define IDL_AX_NONE 4
#define IDL_AX_NOBOX 8
#define IDL_AX_NOZERO 16

typedef struct {		/* System variable for axis */
  IDL_STRING title;		/* Axis title */
  int type;			/* 0 = normal linear, 1=log. */
  int style;			/* 0 = norm, AX_EXTEND, AX_EXACT, AX_NONE,
				   AX_NOBOX. */
  int nticks;			/* # of ticks, 0=auto, -1 = none */
  float ticklen;		/* Tick length, normalized */
  float thick;			/* Axis thickness */
  float range[2];		/* Min and max of endpoints */
  float crange[2];		/* Current min & max */
  float s[2];			/* Scale factors, screen = data*s[1]+s[0] */
  float margin[2];		/* Margin size, in char units. */
  float omargin[2];		/* Outer margin, in char units */
  float window[2];		/* data WINDOW coords, normal units */
  float region[2];		/* Plot region, normal units */
  float charsize;		/* Size of annotations */
  int minor_ticks;		/* Minor ticks */
  float tickv[IDL_MAX_TICKN];	/* Position of ticks */
  IDL_STRING annot[IDL_MAX_TICKN];   /* Annotation */
  IDL_LONG gridstyle;		/* tick linestyle */
  IDL_STRING format;		/* Axis label format/procedure */
  /* After here, the elements are not accessible to the user via
   the system variables: */
  IDL_VPTR ret_values;		/* Returned tick values */
  int log_minor_ticks;		/* true if the minor tickmarks should be
			           laid out via log10() */
} IDL_AXIS;

/* Define cursor function codes: */

#define IDL_CURS_SET 1		/* Set cursor */
#define IDL_CURS_RD  2		/* Read cursor pos */
#define IDL_CURS_RD_WAIT 3	/* Read cursor with wait */
#define IDL_CURS_HIDE 4		/* Disable cursor */
#define IDL_CURS_SHOW 5		/* Display cursor */
#define IDL_CURS_RD_MOVE 6	/* Read & wait for movement or button */
#define IDL_CURS_RD_BUTTON_UP 7	  /* Wait for button up transition */
#define IDL_CURS_RD_BUTTON_DOWN 8   /* Wait for button down transition */
#define IDL_CURS_HIDE_ORIGINAL 9   /* Restore cursor to its original shape
				      (window systems) instead of blanking
				      it. */

/* Define coordinate system types: */
#define IDL_COORD_DATA 0
#define IDL_COORD_DEVICE 1
#define IDL_COORD_NORMAL 2
#define IDL_COORD_MARGIN 3
#define IDL_COORD_IDEVICE 4

#define IDL_PX 0		/* Subscripts for each point member */
#define IDL_PY 1
#define IDL_PZ 2
#define IDL_PH 3

typedef union {			/* Describe a point: */
  struct {			/* Discrete point, ref by .x, .y, or .z: */
    float x,y,z,h;		/* Homogenous coordinates */
  } d;
  struct {			/* Integer discrete: */
    int x,y,z,h;
  } i;
  float p[4];			/* Point refered to by [0],[1],... for
				   x, y, etc.*/
  int ip[4];			/* Integer representation, only valid
				   for COORD_IDEVICE. */
} IDL_GR_PT;

typedef struct {
	IDL_GR_PT origin;
	IDL_GR_PT size;
} IDL_GR_BOX;
    
typedef struct {		/* Attributes structure for points & lines */
  IDL_ULONG color;		/* Specifys all that can go wrong w/ graphic */
  float thick;
  int linestyle;
  float *t;			/* NULL for no 3d transform, or pointer to
				   4 by 4 matrix. */
  int *clip;			/* NULL for no clipping or ^ to [2][2]
				   clipping rectangle in device coord. */
  IDL_AXIS *ax,*ay,*az;		/* Axis definitions */
  int chl;			/* For devices with multiple channels */
} IDL_ATTR_STRUCT;

typedef struct {		/* Graphic text attribute structure.
				   Passed to text routines. */
  int font;			/* 0=hdw, -1=hershey, 1=TrueType */
  int axes;			/* Text axes, 0 = xy, 1 = xz, 2 =yz,
				   3 = yx, 4 = zx, 5 = zy */
  float size;			/* Text size, 1.0 = normal */
  float orien;			/* Orientation, degrees CCW from normal */
  float align;			/* Justification, 0.0 = left, 1.0 = right,
				   0.5  = centered. */
} IDL_TEXT_STRUCT;


/* Structure that defines secondary paramenters for imaging */
typedef struct {		/* Imaging attribute structure */
  short xsize_exp;		/* Non-0 if xsize is EXPlictily set by user */
  short ysize_exp;		/* Non-0 if ysize is EXPlictily set by user */
  IDL_LONG xsize, ysize;	/* Requested size of image (dev coords) */
  int chl;			/* Channel */
  int order;			/* Image order - 0 bottom to top */
  /* Three element array giving the stride between colors of the same
     pixel, adjacent pixels of the same color, and rows of the same color.
     color_stride[0] is non-zero for true color. */
  int color_stride[3];
  int image_is_scratch;		/* True if source image is a temp */
  int b_per_pixel;		/* # of bytes/pixel */
} IDL_TV_STRUCT;



typedef struct {		/* Structure defining polygon fills */
  enum {
    POLY_SOLID, POLY_PATTERN, POLY_IMAGE, POLY_GOURAUD, POLY_IMAGE3D
    } fill_type;
  IDL_ATTR_STRUCT *attr;		/* Graphics attribute structure */
  IDL_PRO_PTR routine;		/* Drawing routine to use */
    union {			/* Operation dependent params: */
    struct {			/* Image fill  */
      UCHAR *data;		/* Fill data for image fill */
      int d1, d2;		/* Dimensions of fill data */
      float *im_verts;		/* Image coords of verts */
      float *im_w;		/* Screen vert W coords */
      UCHAR interp;		/* TRUE to interpolate in image space */
      UCHAR transparent;	/* Transparency threshold, 0 for none */
    } image;
    struct {			/* Line-pattern fill: */
      float angle;		/* Fill orientation in degrees */
      int spacing;		/* Line spacing, in device units */
      float ct, st;		/* Cos / sin of rotation */
    } lines;
    int fill_style;		/* Hardware dependent fill style for
				   POLY_SOLID */
  } extra;
  struct {			/* Info used only with Z buffer device */
    float *z;			/* The Z values */
    int *shades;		/* Shading values at verts for POLY_GOURAUD */
  } three;
} IDL_POLYFILL_ATTR;


typedef struct {		/* Struct containing last mouse status */
  int x,y;			/* X & Y device coordinates */
  int button;			/* Button status bits */
  int time;			/* Time stamp,  not present in all devices */
} IDL_MOUSE_STRUCT;

/*
 * IDL_DEVICE_CORE defines the core functions required by every
 * device driver. Most fields can be filled with
 * a NULL indicating the ability dosen't exist. draw and erase are
 * exceptions to this --- If you can't do that much, why bother with a driver?
*/
typedef struct {
  void (* draw)(IDL_GR_PT *p0, IDL_GR_PT *p1, IDL_ATTR_STRUCT *a);
  int (* text)(IDL_GR_PT *p, IDL_ATTR_STRUCT *ga, IDL_TEXT_STRUCT *ta,
	       char *text);
  void (* erase)(IDL_ATTR_STRUCT *a);	/* erase */
				/* cursor inquire and set */
  void (* cursor)(int funct, IDL_MOUSE_STRUCT *m);
				/* Fill irregular polygon */
  void (* polyfill)(int *x, int *y, int n, IDL_POLYFILL_ATTR *poly);
				/* Returning to interactive mode */
  void (* inter_exit)(void);
  void (* flush)(void);		/* Flush entry */
  void (* load_color)(IDL_LONG start, IDL_LONG n);
				/* Pixel input/output */
  void (* rw_pixels)(UCHAR *data, int x0, int y0, int nx,
		     int ny, int dir, IDL_TV_STRUCT *secondary);
				/* DEVICE procedure */
  void (* dev_specific)(int argc, IDL_VPTR *argv, char *argk);
				/* HELP,/DEVICE */
  void (* dev_help)(int argc, IDL_VPTR *argv);
  void (* load_rtn)(void);	/* Call when driver is loaded */
} IDL_DEVICE_CORE;


/*
 * IDL_DEVICE_WINDOW contains pointers to functions that accomplish
 * window system operations. If the device is a window system,
 * every field in this struct must point at a valid function,
 * they're called without checking.
 */
/*
 * Older mips compilers (Ultrix 4.2) can't handle Ansi function prototypes
 * inside a struct, so we typedef each function and use the defined
 * type inside the struct.
 */
typedef struct {		/* Procedures & functions for image device: */
  void (* window_create)(int argc, IDL_VPTR *argv,char *argk);
  void (* window_delete)(int argc, IDL_VPTR *argv);
  void (* window_show)(int argc, IDL_VPTR *argv, char *argk);
  void (* window_set)(int argc, IDL_VPTR *argv);
  IDL_VPTR (* window_menu)(int argc, IDL_VPTR *argv, char *argk);
} IDL_DEVICE_WINDOW;

/*
 * IDL_DEVICE_DEF is the interface between a device driver and the rest
 * of the Structure defining a graphics device. Every field in this
 * structure must contain valid information --- it is used without
 * any error checking.
 */
typedef struct {		/* Device descriptor, mostly static attributes
				   and definitions: */
  IDL_STRING name;		/* Device name */
  int t_size[2];		/* Total size in device coordinates */
  int v_size[2];		/* Visible area size, device coords */
  int ch_size[2];		/* Default character sizes */
  float px_cm[2];		/* Device units / centimeter, x & y. */
  int n_colors;			/* # of possible simultaneous colors */
  int table_size;		/* # of color table elements */
  int fill_dist;		/* minimum line spacing for solid fill */
  int window;			/* Current window number */
  int unit;			/* Unit number of output file */
  int flags;			/* Advertise limitations and abilities */
  int origin[2];		/* Display XY (pan/scroll) origin */
  int zoom[2];			/* Display XY zoom factors */
  float aspect;			/* Aspect ratio, = v_size[0] / v_size[1]. */
  IDL_DEVICE_CORE core;		/* Core graphics */
  IDL_DEVICE_WINDOW winsys;	/* Window system. Only required if D_WINDOWS */
  char *reserved;		/* Set to zero. */
} IDL_DEVICE_DEF;


/* Define bits in IDL_DEVICE_DEF flags: */

#define IDL_D_SCALABLE_PIXELS 1	  /* True if pixel size is variable (e.g. PS)*/
#define IDL_D_ANGLE_TEXT (1 << 1)   /* True if device can output text at
				       angles */
#define IDL_D_THICK (1 << 2)	/* True if line thickness can be set */
#define IDL_D_IMAGE (1 << 3)	/* True if capable of imaging */
#define IDL_D_COLOR (1 << 4)	/* True if device supports color */
#define IDL_D_POLYFILL (1 << 5)	  /* True if device can do polyfills */
#define IDL_D_MONOSPACE (1<<6)	 /* True if device has only monspaced text */
#define IDL_D_READ_PIXELS (1<<7)   /* True if device can read back pixels */
#define IDL_D_WINDOWS (1<<8)	/* True if device supports windows */
#define IDL_D_WHITE_BACKGROUND (1<<9)	/* True if device background is
					   white, like PostScript. */
#define IDL_D_NO_HDW_TEXT (1<<10)   /* True if device has no hardware text */
#define IDL_D_POLYFILL_LINE (1<<11)   /* True to use device driver for line
					 style polyfills. */
#define IDL_D_HERSH_CONTROL (1<<12)   /* True if device accepts hershey style
					 control characters. */
#define IDL_D_PLOTTER (1<<13)	/* True if pen plotter */
#define IDL_D_WORDS (1<<14)	/* True if device images can be words */
#define IDL_D_KANJI (1 << 15)	/* Device has Kanji characters */
#define IDL_D_WIDGETS (1 << 16)	  /* Device supports graphical user
				     interfaces */
#define IDL_D_Z (1 << 17)	/* Device is 3d */
#define IDL_D_TRUETYPE_FONT (1 << 18) /* Device supports TrueType fonts. */

typedef struct {
  int background;		/* Background color */
  float charsize;		/* Global Character size */
  float charthick;		/* Character thickness */
  int clip[6];			/* Clipping rectangle, normalized coords */
  IDL_ULONG color;		/* Current color */
  int font;			/* Font */
  int linestyle;		/* Line style */
  int  multi[5];		/* Cnt, Cols/rows, major dir for multi plts. */
  int clip_off;			/* True if clipping is disabled */
  int noerase;			/* No erase flag */
  int nsum;			/* Number of points to sum */
  float position[4];		/* Default window */
  int psym;			/* Marker symbol */
  float region[4];		/* Default plotting region */
  IDL_STRING subtitle;		/* Plot subtitle */
  float symsize;		/* Symbol size */
  float t[16];			/* Matrix (4x4) of homogeneous transform */
  int  t3d_on;			/* True if 3d homo transform is on */
  float thick;			/* Line thickness */
  IDL_STRING title;		/* Main plot title */
  float ticklen;		/* Tick length */
  int chl;			/* Default channel */
  IDL_DEVICE_DEF *dev;		/* Current output device, not user accesible */
} IDL_PLOT_COM;



#endif				/* graphics_IDL_DEF */




/***** Definitions from keyword *****/


#ifndef keyword_IDL_DEF
#define keyword_IDL_DEF

/* Bit values of IDL_KW_PAR flags field: */

#define IDL_KW_ARRAY (1 << 12)
/* If specified array is required, otherwise scalar required */

#define IDL_KW_OUT (1 << 13)
/* Indicates parameter is an OUTPUT parameter passed by reference.
   Expressions are excluded.  The address of the IDL_VARIABLE is stored in
   the value field.  Otherwise, no checking is performed.  Special hint:
   to find out if a IDL_KW_OUT parameter is specified, use 0 for the type,
   and IDL_KW_OUT | IDL_KW_ZERO for the flags.  The value field will either
   contain NULL or the pointer to the variable. */


#define IDL_KW_VIN (IDL_KW_OUT | IDL_KW_ARRAY)
/* Parameter is an INPUT parameter passed by reference.  The address
   of the IDL_VARIABLE or expression is stored in the value field as with
   IDL_KW_OUT.  If this flag is specified, kw_cleanup() must be called to
   properly reap temporaries that may have been allocated. */


#define IDL_KW_ZERO (1 << 14)
/* If set, zero the parameter before parsing the keywords.  I.e. if
   this bit is set, and the parameter is not specified, the value will
   always be 0. */

#define IDL_KW_VALUE (1 << 15)
/* If this bit is set and the keyword is present, and its value is
   non-zero, the low 12 bits of this field will be inclusive 'or'ed with
   the longword pointed to by IDL_KW_PAR.value.  Be sure that the type field
   contains TYP_LONG.  The largest value that may be specified is
   (2^12)-1.  Negative values are not allowed.  For example, if the
   IDL_KW_PAR struct contains:

   "DEVICE", TYP_LONG, 1, IDL_KW_VALUE | 4 |IDL_KW_ZERO,NULL,&(char *) xxx,
   "NORMAL", TYP_LONG, 1, IDL_KW_VALUE | 3, NULL, &(char *) xxx,

   then xxx will contain a 3 if /NORMAL, or NORMAL = (expr) is
   present, a 4 if /DEVICE is present, 7 if both are set, and 0 if
   neither.  IDL_KW_ZERO can also be used in combination with this flag, use
   it only once for each IDL_KW_PAR.value.  */



#define IDL_KW_VALUE_MASK ((1 << 12) -1)   /* Mask for value part */

/* Use IDL_KW_FAST_SCAN as the first element of the keyword array if
   there are more than approximately 5 or 10 elements in the keyword
   array.  The IDL_KW_PAR structure defined by this macro is used to
   point to a list of elements to zero, and speeds processing of long
   keyword lists.  NEVER touch the contents of this structure.
   */
#define IDL_KW_FAST_SCAN { (char *)"", 0,0,0,0,0 }


typedef struct {
  char *keyword;		/* ^ to Keyword string, NULL terminated.
				   A NULL keyword string pointer value
				   terminates the keyword structure.  Strings
				   must be UPPER case and in LEXICAL order .*/
  UCHAR type;			/* Type of data required.  For scalars the
				   only allowable types are TYP_STRING,
				   TYP_LONG, TYP_FLOAT and
				   TYP_DOUBLE. For arrays, this may be
				   any simple type or 0 for no conversion. */
  unsigned short mask;		/* Enable mask.  This field is AND'ed with
				   the mask field in the call to
				   GET_IDL_KW_PARAMS, and if the result is
				   non-zero the keyword is used. If it is 0,
				   the keyword is ignored.  */
  unsigned short flags;		/* Contains  flags as described above */
  int *specified;		/* Address of int to set on return if
				   param is specified.  May be null if this
				   information is not required. */
  char *value;			/* Address of value to return.  In the
				   case of arrays,  this value points to
				   the IDL_KW_ARR_DESC structure for the data
				   to be returned.
				   */
} IDL_KW_PAR;

typedef struct {		/* Descriptor for array's that are returned */
  char *data;			/* Address of array to receive data. */
  IDL_MEMINT nmin;		/* Minimum # of elements allowed. */
  IDL_MEMINT nmax;		/* Maximum # of elements allowed. */
  IDL_MEMINT n;			/* # present, (Returned value). */
} IDL_KW_ARR_DESC;



#define IDL_KW_CLEAN_ALL 0	/* Codes for kw_cleanup, clean all temps
				   and strings */
#define IDL_KW_MARK 1		/* Mark string stack before calling get_
				   kw_params.  */
#define IDL_KW_CLEAN 2		/* Clean temps & strings created since
				   last call with KW_MARK. */

#endif				/* keyword_IDL_DEF */




/***** Definitions from lmgr *****/

/*
 * Define licensing request codes for the exported API
 */
#define IDL_LMGR_CLIENTSERVER     0x01
#define IDL_LMGR_DEMO             0x02
#define IDL_LMGR_EMBEDDED         0x04
#define IDL_LMGR_RUNTIME          0x08
#define IDL_LMGR_STUDENT          0x10
#define IDL_LMGR_TRIAL            0x20
#define IDL_LMGR_CALLAPPNOCHECKOUT      0x40
#define IDL_LMGR_CALLAPPLICINTERNAL		0x80




/***** Definitions from os *****/

#ifndef os_IDL_DEF
#define os_IDL_DEF

/* Structure passed to IDL_GetUserInfo() */
typedef struct {
     char *logname;			/* Users login name */
     char host[64];			/* The machine name */
     char wd[IDL_MAXPATH+1];		/* The current directory */
     char date[25];			/* The current date */
   } IDL_USER_INFO;

/* SPAWN allows more arguments on the Macintosh than on the other platforms. */
#ifdef MAC
#define IDL_MAXSPAWNPARAMS IDL_MAXPARAMS
#else
#define IDL_MAXSPAWNPARAMS 2
#endif

#ifdef VMS
#include <descrip.h>
#endif

#endif				/* os_IDL_DEF */




/***** Definitions from pout *****/

#ifndef pout_IDL_DEF
#define pout_IDL_DEF

/*** Mask values for flags argument to pout() ***/
#define IDL_POUT_SL         1	/* Start on a new line */
#define IDL_POUT_FL         2	/* Finish current line */
#define IDL_POUT_NOSP       4	/* Don't add leading space */
#define IDL_POUT_NOBREAK    8	/* Don't start a new line if too long */
#define IDL_POUT_LEADING    16	/* Print leading text at start of line */
#define IDL_POUT_GET_POS    32	/* Get current file position */
#define IDL_POUT_SET_POS    64	/* Set current file position */


/*** Structure for control argument to pout() ***/
typedef struct {
  int unit;			/* LUN of open file */
  int curcol;			/* Current output column */
  int wrap;			/* # chars at which buf should flush */
  char *leading;		/* String to output at start of each line */
  int leading_len;		/* Length of leading w/o terminating null */
  char *buf;			/* ^ to output buffer. Must be max_len chars */
  int max_len;			/* Length of buffer */
} IDL_POUT_CNTRL;

#endif				/* pout_IDL_DEF */




/***** Definitions from raster *****/

#ifndef raster_IDL_DEF
#define raster_IDL_DEF


/*** Allowed values for dither_method field of IDL_RASTER_DEF struct ***/
#define IDL_DITHER_REVERSE          0 /* "Un-dither" back to bytes */
#define IDL_DITHER_THRESHOLD        1 /* Threshold dithering  */
#define IDL_DITHER_FLOYD_STEINBERG  2 /* Floyd Steinberg method */
#define IDL_DITHER_ORDERED          3 /* Ordered dither */

/* Values for flags field: */
#define IDL_DITHER_F_WHITE 1	/* Device has white background, dithering
				   module then sets the black bits. */

/*** Convenience values for the bit_tab array of IDL_RASTER_DEF struct ***/
#define IDL_RASTER_MSB_LEFT { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 }
#define IDL_RASTER_MSB_RIGHT { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 }

typedef struct {		/* Information that characterizes a raster */
  UCHAR *fb;			/* Address of frame buffer */
  int nx, ny;			/* Size of frame buffer in pixels */
  int bytes_line;		/* # of bytes per scan line, must be whole. */
  int byte_padding;		/* Pad lines to a multiple of this amount,
				   Must be a power of 2.  */
  int dot_width;		/* The length of a dot for the linestyles.
				   Default = 1. */
  int dither_method;		/* Dithering method code. */
  int dither_threshold;		/* Threshold value for threshold dither */
  UCHAR bit_tab[8];		/* Table of set bits, bit_tab[0] is leftmost,
				   bit_tab[7] is right most. */
  int flags;			/* Raster flags, see above */
} IDL_RASTER_DEF;

#endif				/* raster_IDL_DEF */




/***** Definitions from rline *****/

#ifndef rline_IDL_DEF
#define rline_IDL_DEF

/**** Flags to OR together for options parameter to IDL_Rline() ****/
#define IDL_RLINE_OPT_NOSAVE        1   /* Don't save in recall buffer */
#define IDL_RLINE_OPT_NOJOURNAL     2   /* Don't journal */
#define IDL_RLINE_OPT_JOURCMT       4   /* Put a '; ' at start in journal */
#define IDL_RLINE_OPT_NOEDIT        8   /* Like (!EDIT_INPUT = 0) for one call*/

#endif				/* rline_IDL_DEF */




/***** Definitions from sig *****/

#ifndef sig_IDL_DEF
#define sig_IDL_DEF


#include <signal.h>

/* Dialect confusion from HP-UX */
#if defined(SIGWINDOW) && !defined(SIGWINCH)
#define SIGWINCH SIGWINDOW
#endif

/*
 * Signal sets are represented by this opaque type. The type and length
 * have been selected to be suitable for any platform.
 */
typedef struct {
#ifdef linux
  unsigned long set[_SIGSET_NWORDS];
#else
  double set[4];
#endif
} IDL_SignalSet_t;

/* The IDL definition for all signal handler functions. */
typedef void (* IDL_SignalHandler_t) IDL_ARG_PROTO((int signo));

#endif				/* sig_IDL_DEF */




/***** Definitions from structs *****/

#ifndef structs_IDL_DEF
#define structs_IDL_DEF


/* Valid bits for flags field of IDL_STRUCT_TAG_DEF */
#define IDL_STD_INHERIT	1	/* Type must be a structure. This flag
				   indicates that the structure is inherited
				   (inlined) instead of making it a
				   sub-structure as usual. */

typedef struct {		/* A tag definition for K_MakeSTruct */
  char *name;			/* Name of the tag. Must be upper case
				   and obey the rules for IDL identifiers.
				   In the case of inherited structures, this
				   can be NULL if type is set. Otherwise, it
				   is the name of the structure being inherited
				   and IDL will call a __DEFINE procedure
				   based on that name to define it. */
  IDL_MEMINT *dims;		/* NULL pointer for a scalar tag, otherwise
				   an array giving the array dimensions. The
				   first element is the number of dimensions,
				   and is followed in order by the dimensions
				   themselves. */
  void *type;			/* This may be either a pointer to another
				   structure definition, or a simple IDL
				   type code (IDL_TYP_*) cast to void
				   (e.g. (void *) IDL_TYP_BYTE). If this
				   field is NULL, it indicates that IDL
				   should search for a structure of the
				   given name and fill in the pointer to
				   its structure definition. */
  UCHAR flags;			/* Bitmask made up of IDL_STD_* bits */
} IDL_STRUCT_TAG_DEF;
#endif				/* structs_IDL_DEF */




/***** Definitions from sysnames *****/

#ifndef sysnames_IDL_DEF
#define sysnames_IDL_DEF

/* Structure used for IDL_SysvVersion global variable */
typedef struct {
  IDL_STRING arch;		/* Machine architecture */
  IDL_STRING os;		/* Operating System */
  IDL_STRING os_family;		/* Operating System family
				   (e.g. Unix vs SunOS) */
  IDL_STRING release;		/* Software release */
  IDL_STRING build_date;	/* Date on which this executable was built */
} IDL_SYS_VERSION;

/* Structure used for IDL_SysvErrorState global variable */
typedef struct {
  IDL_STRING name;		/* Symbolic name of current error */
  IDL_STRING block;		/* Name of error block for current error */
  /* Error code. !ERR is often changed by various parts of IDL for
   * non-error reasons. These unfortunate side effects are historical in
   * nature, and cannot be eliminated. !ERROR_STATE.CODE is just like !ERR
   * except that its purpose is pure, it only contains the code of the
   * last error, and always matches the rest of !ERROR_STATE.
   */
  IDL_LONG code;		
  IDL_LONG sys_code[2];		/* System error code */
  IDL_STRING msg;		/* Text of IDL error message */
  IDL_STRING sys_msg;		/* System component of error message */
  IDL_STRING msg_prefix;	/* Prefix attached to all error messages */
} IDL_SYS_ERROR_STATE;

/*
 * These #defines allow use of older error related system variables and
 * map them to the corresponding field in !ERROR_STATE. These #defines
 * will eventually be obsoleted and moved to obsolete.h. Programmers
 * should convert their code to use the new name in preparation for this.
 */
#define IDL_SysvErrString IDL_SysvErrorState.msg
#define IDL_SysvSyserrString IDL_SysvErrorState.sys_msg
#define IDL_SysvErrorCode IDL_SysvErrorState.code
#define IDL_SysvSyserrorCodes IDL_SysvErrorState.sys_code

#endif				/* sysnames_IDL_DEF */




/***** Definitions from tout *****/

#ifndef tout_IDL_DEF
#define tout_IDL_DEF

typedef void (* IDL_TOUT_OUTF)(int flags, char *buf, int n);

#define IDL_TOUT_F_STDERR   1	/* Output to stderr instead of stdout */
#define IDL_TOUT_F_NLPOST   4	/* Output a newline at end of line */
     
#endif				/* tout_IDL_DEF */




/***** Definitions from uicb *****/

typedef int (* IDL_UicbMacEvent_t)(void *event);





/***** Definitions from ur_main *****/

#ifndef ur_main_IDL_DEF
#define ur_main_IDL_DEF


/* Values that are OR'd together to form the options argument to IDL_Init() */
#define IDL_INIT_GUI	1	/* Use the GUI interface. */
#define IDL_INIT_GUI_AUTO (IDL_INIT_GUI|2)
				/* Try to use a GUI if possible. If that
				   fails and the OS supports ttys, use
				   the standard tty interface. Note that
				   this code includes IDL_INIT_GUI. */
#define IDL_INIT_RUNTIME 4	/* RunTime IDL. */
#define IDL_INIT_EMBEDDED (IDL_INIT_RUNTIME|8)
				/* Embedded IDL. Note that this code includes
				   IDL_INIT_RUNTIME. */
#define IDL_INIT_NOLICALIAS 16	 /* Our FlexLM (Unix/VMS) floating license
				    policy is to alias all IDL sessions that
				    share the same user/system/display to the
				    same license. If no_lic_alias is set,
				    this IDL session will force a unique
				    license to be checked out. In this case,
				    we allow the user to change the DISPLAY
				    environment variable. This is useful for
				    RPC servers that don't know where their
				    output will need to go before invocation.*/
#define IDL_INIT_BACKGROUND 32	/* This tells IDL that it is going to be used
				   in a background mode by some other program,
				   and it is not in control of the user's
				   input command line. One effect of this is
				   that XMANAGER will block, since IDL cannot
				   use its active command line functionality
				   to dispatch events.

				   Normally under Unix, if IDL sees that
				   stdin and stdout are ttys, it puts the tty
				   into raw mode and uses termcap/terminfo to
				   handle command line editing. When using
				   callable IDL in a background process that
				   isn't doing I/O to the tty, the termcap
				   initialization can cause the process
				   to block because of job control from the
				   shell with a message like "Stopped (tty
				   output) idl". Setting this option prevents
				   all tty edit functions and disables the
				   calls to termcap. I/O to the tty is done
				   with a simple fgets()/printf().
				   In the case of IDL_INIT_GUI, this is
				   ignored. */

#define IDL_INIT_NOTTYEDIT IDL_INIT_BACKGROUND
				/* Renamed it to better reflect it's more
				   general functionality. */

#define IDL_INIT_QUIET 64	/* Suppresses the startup announcement and
				   message of the day. */

#define IDL_INIT_STUDENT 128	/* IDL Student Edition */

#define IDL_INIT_CLIENT  256
#define IDL_INIT_GUI_CLIENT (IDL_INIT_CLIENT|IDL_INIT_GUI)
                                /* Start IDLDE if it wasn't started, else
				   just send the filename specified to
				   already running IDLDE. */

#define IDL_INIT_DEMO	512	/* Force IDL into demo mode */

#define IDL_INIT_VAX_FLOAT 1024   /* VMS-only: Cause the /VAX_FLOAT
				     keyword to OPEN and CALL_EXTERNAL
				     default to TRUE instead of the usual
				     FALSE. This lets VAX dependant code
				     run as if VMS/IDL still used VAX
				     floating point. */

#endif				/* main_IDL_DEF */




/***** Definitions from widgets *****/

#ifndef widgets_IDL_DEF
#define widgets_IDL_DEF

typedef void (* IDL_WIDGET_STUB_SET_SIZE_FUNC)
     IDL_ARG_PROTO((IDL_ULONG id, int width, int height));

#endif  /* widgets_DEF */




/***** Definitions from zfiles *****/

#ifndef zfiles_IDL_DEF
#define zfiles_IDL_DEF
/**** Access field bits in IDL_FILE_DESC and IDL_FILE_STAT ****/
#define IDL_OPEN_R          1	/* Open file for reading */
#define IDL_OPEN_W          2	/* Open file for writing */
#define IDL_OPEN_NEW        4	/* Unix - Truncate old file contents.
                                   VMS - Use a new file. */
#define IDL_OPEN_APND       8	/* File open with pointer at EOF */

/**** Flags field bits in IDL_FILE_DESC and IDL_FILE_STAT ****/
#define IDL_F_ISATTY        1	/* Is a terminal */
#define IDL_F_ISAGUI	    2	/* Is a Graphical User Interface */
#define IDL_F_NOCLOSE       4	/* Don't let user close */
#define IDL_F_MORE          8	/* Use more(1) like pager for fmt output */
#define IDL_F_XDR           16	 /* Is a XDR file */
#define IDL_F_DEL_ON_CLOSE  32	 /* Delete on close */
#define IDL_F_SR            64	 /* Is a SAVE/RESTORE file. */
#define IDL_F_SWAP_ENDIAN   128	/* File has opposite byte order than current
				   system. */
#define IDL_F_VAX_FLOAT	    (1 << 8) /* Binary float and double are in VAX
					F and D format. Implies that file is
					little endian. */
#define IDL_F_UNIX_F77      (1 << 9) /* Unformatted f77(1) I/O */
#define IDL_F_UNIX_PIPE     (1 << 10) /* fptr is to a socketpair(2) */
#define IDL_F_UNIX_NOSTDIO (1 << 11) /* Call read(2) and write(2) directly */
#define IDL_F_UNIX_SPECIAL  (1 << 12) /* It's a device/special file */
#define IDL_F_VMS_FIXED     (1 << 13) /* Fixed length records */
#define IDL_F_VMS_VARIABLE  (1 << 14) /* Variable length records */
#define IDL_F_VMS_SEGMENTED (1 << 15) /* FORTRAN segmented var len records */
#define IDL_F_VMS_STREAM    (1 << 16) /* Stream file */
				/* When reading a non-stream file via
				   VMS stdio, there are two possible
				   approaches. One is to simply read
				   the file as a stream and let the RMS
				   stuff show. The other is to try to
				   re-write the data into a "logical
				   data stream". Normally, IDL takes
				   the second approach because it
				   allows  code to work easily between Unix and
				   VMS. For the user OPEN though, the
				   first approach is better because it
				   is more robust and avoids RMS
				   buffer size limitations. the
				   STREAM_STRICT modifier flag is used
				   in this case. */
#define IDL_F_VMS_STREAM_STRICT (1 << 17)
#define IDL_F_VMS_RMSBLK    (1 << 18)   /* RMS Block Mode access */
#define IDL_F_VMS_RMSBLKUDF (1 << 19)   /* RMS block mode files are created
					   with FIXED length 512 byte records.
					   This bit indicates that the RMS
					   UNDEFINED record type should be
					   used. */
#define IDL_F_VMS_INDEXED    (1 << 20)   /* Indexed file */
#define IDL_F_VMS_PRINT      (1 << 21)   /* Send to SYS$PRINT on close */
#define IDL_F_VMS_SUBMIT     (1 << 22)   /* Send to SYS$BATCH on close */
#define IDL_F_VMS_TRCLOSE    (1 << 23)   /* Truncate file allocation on close*/
#define IDL_F_VMS_CCLIST     (1 << 24)   /* CR/LF carriage control */
#define IDL_F_VMS_CCFORTRAN  (1 << 25)   /* FORTRAN style carriage control */
#define IDL_F_VMS_CCNONE     (1 << 26)   /* Explicit carriage control */
#define IDL_F_VMS_SHARED     (1 << 27)   /* Shared access */
#define IDL_F_VMS_SUPERSEDE  (1 << 28)   /*Supersede existing version on open*/
#define IDL_F_DOS_NOAUTOMODE (1 << 29)   /* Don't switch the mode */
#define IDL_F_DOS_BINARY     (1 << 30)   /* File is in binary mode (^J) */
  
/* Sets the IDL_F_NOCLOSE bit for file unit. */
#define IDL_FILE_NOCLOSE(unit) IDL_FileSetClose((unit), FALSE)

/* Clear the IDL_F_NOCLOSE bit for file unit. */
#define IDL_FILE_CLOSE(unit) IDL_FileSetClose((unit), TRUE)

/**** File units that map to standard units ****/
#define IDL_STDIN_UNIT      0
#define IDL_STDOUT_UNIT     -1
#define IDL_STDERR_UNIT     -2
#define IDL_NON_UNIT        -100    /* Gauranteed to be an invalid unit */

/* Valid flags to bit-OR together for IDL_FileEnsureStatus() flags argument */
#define IDL_EFS_USER        1       /* Must be user unit (1 - MAX_USER_FILES) */
#define IDL_EFS_OPEN        2       /* Unit must be open */
#define IDL_EFS_CLOSED      4       /* Unit must be closed */
#define IDL_EFS_READ        8       /* Unit must be open for input */
#define IDL_EFS_WRITE       16      /* Unit must be open for output */
#define IDL_EFS_NOTTY       32      /* Unit cannot be a tty */
#define IDL_EFS_NOGUI       64      /* Unit cannot be a tty */
#define IDL_EFS_NOPIPE      128      /* Unit cannot be a pipe */
#define IDL_EFS_NOXDR       256     /* Unit cannot be a XDR file */
#define IDL_EFS_ASSOC       512     /* Unit can be assoc'd. This implies USER,
                                   OPEN, NOTTY, NOPIPE, and NOXDR, in
                                   addition to other OS specific concerns */
#define IDL_EFS_NOT_NOSTDIO 1024	/* Under Unix, file wasn't opened with
				   IDL_F_UNIX_NOSTDIO attribute. */


/**** Struct for global variable term, filled by IDL_FileInit() ****/
typedef struct {
#ifdef IDL_OS_HAS_TTYS
  char *name;                   /* Name of terminal type */
  char is_tty;                  /* True if stdin is a terminal */
#endif
  int lines;                    /* Lines on screen */
  int columns;                  /* Width of output */
} IDL_TERMINFO;


/**** Struct that is filled in by IDL_FileStat() ****/
typedef struct {
  char *name;
  short access;
  IDL_LONG flags;
  FILE *fptr;
  struct {
    unsigned short mrs;
  } rms;
} IDL_FILE_STAT;

#endif				/* zfiles_IDL_DEF */




/***** Definitions from ztimer *****/

#ifndef ztimer_IDL_DEF
#define ztimer_IDL_DEF

typedef void (* IDL_TIMER_CB)();
#ifdef VMS
typedef long IDL_TIMER_CONTEXT;
#else
typedef void (* IDL_TIMER_CONTEXT)();
#endif
typedef IDL_TIMER_CONTEXT *IDL_TIMER_CONTEXT_PTR;

#endif				/* ztimer_IDL_DEF */




/* Forward declarations for all exported routines and data */


extern char *IDL_CDECL IDL_FilePathFromRoot IDL_ARG_PROTO((char *pathbuf,
        IDL_STRING *root, char *file, char *ext, int nsubdir,  char
        **subdir));
extern char *IDL_CDECL IDL_FilePath IDL_ARG_PROTO((char *pathbuf, char
        *file, char *ext, int nsubdir, char **subdir));
extern void *IDL_CDECL IDL_MemAlloc IDL_ARG_PROTO((IDL_MEMINT n, char
        *err_str, int action));
extern void IDL_CDECL IDL_MemFree IDL_ARG_PROTO((IDL_REGISTER void *m,
        char *err_str, int action));
extern void *IDL_CDECL IDL_MemAllocPerm IDL_ARG_PROTO((IDL_MEMINT n, char
        *err_str, int action));
extern void IDL_CDECL IDL_TimerSet IDL_ARG_PROTO((IDL_LONG length,
        IDL_TIMER_CB callback, int from_callback, IDL_TIMER_CONTEXT_PTR
        context));
extern void IDL_CDECL IDL_TimerCancel IDL_ARG_PROTO((IDL_TIMER_CONTEXT
        context));
extern void IDL_CDECL IDL_TimerBlock IDL_ARG_PROTO((int stop));
extern IDL_UicbMacEvent_t IDL_CDECL IDL_UicbRegMacEvent
        IDL_ARG_PROTO((IDL_UicbMacEvent_t func));
extern void IDL_CDECL IDL_SignalSetInit IDL_ARG_PROTO((IDL_SignalSet_t
        *set, int signo));
extern void IDL_CDECL IDL_SignalSetAdd IDL_ARG_PROTO((IDL_SignalSet_t
        *set, int signo));
extern void IDL_CDECL IDL_SignalSetDel IDL_ARG_PROTO((IDL_SignalSet_t
        *set, int signo));
extern int IDL_CDECL IDL_SignalSetIsMember IDL_ARG_PROTO((IDL_SignalSet_t
        *set, int signo));
extern void IDL_CDECL IDL_SignalMaskGet IDL_ARG_PROTO((IDL_SignalSet_t
        *set));
extern void IDL_CDECL IDL_SignalMaskSet IDL_ARG_PROTO((IDL_SignalSet_t
        *set, IDL_SignalSet_t *oset));
extern void IDL_CDECL IDL_SignalMaskBlock IDL_ARG_PROTO((IDL_SignalSet_t
        *set, IDL_SignalSet_t *oset));
extern void IDL_CDECL IDL_SignalBlock IDL_ARG_PROTO((int signo,
        IDL_SignalSet_t *oset));
extern void IDL_CDECL IDL_SignalSuspend IDL_ARG_PROTO((IDL_SignalSet_t
        *set));
extern int IDL_CDECL IDL_SignalRegister IDL_ARG_PROTO((int signo,
        IDL_SignalHandler_t func, int msg_action));
extern int IDL_CDECL IDL_SignalUnregister IDL_ARG_PROTO((int signo,
        IDL_SignalHandler_t func, int msg_action));
extern int IDL_CDECL IDL_GetKbrd IDL_ARG_PROTO((int should_wait));
extern void IDL_CDECL IDL_TerminalRaw IDL_ARG_PROTO((int to_from, int
        fnin));
extern void IDL_CDECL IDL_Pout IDL_ARG_PROTO((IDL_POUT_CNTRL *control,
        ...));
extern void IDL_CDECL IDL_PoutRaw IDL_ARG_PROTO((int unit, char *buf, int
        n));
extern IDL_PLOT_COM IDL_PlotCom;
extern UCHAR IDL_ColorMap[];
extern IDL_PLOT_COM *IDL_CDECL IDL_PlotComAddr IDL_ARG_PROTO((void));
extern UCHAR *IDL_CDECL IDL_ColorMapAddr IDL_ARG_PROTO((void));
extern void IDL_CDECL IDL_PolyfillSoftware IDL_ARG_PROTO((int *x, int *y,
        int n, IDL_POLYFILL_ATTR *s));
extern double IDL_CDECL IDL_GraphText IDL_ARG_PROTO((IDL_GR_PT *p,
        IDL_ATTR_STRUCT *ga, IDL_TEXT_STRUCT *a, char *text));
extern void IDL_CDECL IDL_StrDup IDL_ARG_PROTO((IDL_REGISTER IDL_STRING
        *str, IDL_REGISTER IDL_MEMINT n));
extern void IDL_CDECL IDL_StrDelete IDL_ARG_PROTO((IDL_STRING *str,
        IDL_MEMINT n));
extern void IDL_CDECL IDL_StrStore IDL_ARG_PROTO((IDL_STRING *s, char
        *fs));
extern void IDL_CDECL IDL_StrEnsureLength IDL_ARG_PROTO((IDL_STRING *s,
        int n));
extern IDL_VPTR IDL_CDECL IDL_StrToSTRING IDL_ARG_PROTO((char *s));
extern void IDL_CDECL IDL_TerminalRaw IDL_ARG_PROTO((int to_from, int
        fnin));
extern int IDL_CDECL IDL_GetKbrd IDL_ARG_PROTO((int should_wait));
extern int IDL_STDCALL IDL_InitOCX IDL_ARG_PROTO((void *pInit));
extern int IDL_CDECL IDL_GetKbrd IDL_ARG_PROTO((int should_wait));
extern void IDL_CDECL IDL_ExitRegister IDL_ARG_PROTO((IDL_PRO_PTR proc));
extern void IDL_CDECL IDL_WidgetIssueStubEvent IDL_ARG_PROTO((char *rec,
        IDL_LONG value));
extern void IDL_CDECL IDL_WidgetSetStubIds IDL_ARG_PROTO((char *rec,
        unsigned long t_id, unsigned long b_id));
extern void IDL_CDECL IDL_WidgetGetStubIds IDL_ARG_PROTO((char *rec,
        unsigned long *t_id, unsigned long *b_id));
extern void IDL_CDECL IDL_WidgetStubLock IDL_ARG_PROTO((int set));
extern void *IDL_CDECL IDL_WidgetStubGetParent IDL_ARG_PROTO((IDL_ULONG
        id, char *szDisplay));
extern char *IDL_CDECL IDL_WidgetStubLookup IDL_ARG_PROTO((IDL_ULONG id));
extern void IDL_CDECL IDL_WidgetStubSetSizeFunc IDL_ARG_PROTO((char *rec,
        IDL_WIDGET_STUB_SET_SIZE_FUNC func));
extern void IDL_CDECL IDL_CvtVAXToFloat IDL_ARG_PROTO((float *fp,
        IDL_MEMINT n));
extern void IDL_CDECL IDL_CvtFloatToVAX IDL_ARG_PROTO((float *fp,
        IDL_MEMINT n));
extern void IDL_CDECL IDL_CvtVAXToDouble IDL_ARG_PROTO((double *dp,
        IDL_MEMINT n));
extern void IDL_CDECL IDL_CvtDoubleToVAX IDL_ARG_PROTO((double *dp,
        IDL_MEMINT n));
extern void IDL_CDECL IDL_EzCall IDL_ARG_PROTO((int argc, IDL_VPTR argv[],
        IDL_EZ_ARG arg_struct[]));
extern void IDL_CDECL IDL_EzCallCleanup IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[], IDL_EZ_ARG arg_struct[]));
extern void IDL_CDECL IDL_EzReplaceWithTranspose IDL_ARG_PROTO((IDL_VPTR
        *v, IDL_VPTR orig));
extern IDL_MSG_BLOCK IDL_CDECL IDL_MessageDefineBlock IDL_ARG_PROTO((char
        *block_name, int n, IDL_MSG_DEF *defs));
extern void IDL_CDECL IDL_MessageErrno IDL_ARG_PROTO((int code, ...));
extern void IDL_CDECL IDL_MessageErrnoFromBlock
        IDL_ARG_PROTO((IDL_MSG_BLOCK block, int code, ...));
extern void IDL_CDECL IDL_Message IDL_ARG_PROTO((int code, int action,
        ...));
extern void IDL_CDECL IDL_MessageFromBlock IDL_ARG_PROTO((IDL_MSG_BLOCK
        block, int code, int action,...));
extern void IDL_CDECL IDL_MessageVarError IDL_ARG_PROTO((int code,
        IDL_VPTR var, int action));
extern void IDL_CDECL IDL_MessageVarErrorFromBlock
        IDL_ARG_PROTO((IDL_MSG_BLOCK block, int code, IDL_VPTR var, int
        action));
extern void IDL_CDECL IDL_MessageVMS IDL_ARG_PROTO((int code,...));
extern void IDL_CDECL IDL_MessageVMSFromBlock IDL_ARG_PROTO((IDL_MSG_BLOCK
        block, int code,...));
extern int IDL_CDECL IDL_MessageNameToCode IDL_ARG_PROTO((IDL_MSG_BLOCK
        block, char *name));
extern void IDL_CDECL IDL_MessageResetSysvErrorState IDL_ARG_PROTO((void));
extern void IDL_CDECL IDL_MessageSJE IDL_ARG_PROTO((void *value));
extern void *IDL_CDECL IDL_MessageGJE IDL_ARG_PROTO((void));
extern void IDL_CDECL IDL_MessageVE_UNDEFVAR IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOTARRAY IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOTSCALAR IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOEXPR IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOCONST IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOFILE IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOCOMPLEX IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOSTRING IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOSTRUCT IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_REQSTR IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_NOSCALAR IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_STRUC_REQ IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_REQPTR IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_MessageVE_REQOBJREF IDL_ARG_PROTO((IDL_VPTR var,
        int action));
extern void IDL_CDECL IDL_Message_BADARRDNUM IDL_ARG_PROTO((int action));
extern char *IDL_CDECL IDL_Rline IDL_ARG_PROTO((char *s, int n, int unit,
        FILE *stream, int is_tty, char *prompt, int opt));
extern void IDL_CDECL IDL_RlineSetStdinOptions IDL_ARG_PROTO((int opt));
extern void IDL_CDECL IDL_Logit IDL_ARG_PROTO((char *s));
extern int IDL_CDECL IDL_Win32Init IDL_ARG_PROTO((int iOpts, void
        *hinstExe, void *hwndExe, void *hAccel));
extern void IDL_CDECL IDL_ToutPush IDL_ARG_PROTO((IDL_TOUT_OUTF outf));
extern IDL_TOUT_OUTF IDL_CDECL IDL_ToutPop IDL_ARG_PROTO((void));
extern char *IDL_CDECL IDL_VarName IDL_ARG_PROTO((IDL_VPTR v));
extern IDL_VPTR IDL_CDECL IDL_GetVarAddr1 IDL_ARG_PROTO((char *name, int
        ienter));
extern IDL_VPTR IDL_CDECL IDL_GetVarAddr IDL_ARG_PROTO((char *name));
extern IDL_VPTR IDL_CDECL IDL_FindNamedVariable IDL_ARG_PROTO((char *name,
        int ienter));
extern char *IDL_CDECL IDL_MakeTempArray IDL_ARG_PROTO((int type, int
        n_dim, IDL_MEMINT dim[], int init, IDL_VPTR *var));
extern char *IDL_CDECL IDL_MakeTempVector IDL_ARG_PROTO((int type,
        IDL_MEMINT dim, int init, IDL_VPTR *var));
extern void IDL_CDECL IDL_Wait IDL_ARG_PROTO((int argc, IDL_VPTR argv[]));
extern void IDL_CDECL IDL_GetUserInfo IDL_ARG_PROTO((IDL_USER_INFO
        *user_info));
extern void IDL_CDECL IDL_TTYReset IDL_ARG_PROTO((void));
extern short IDL_TapeChl[];
extern short *IDL_CDECL IDL_TapeChlAddr IDL_ARG_PROTO((void));
extern int IDL_CDECL IDL_AddToQueue IDL_ARG_PROTO((char* pString));
extern int IDL_CDECL IDL_GetWait IDL_ARG_PROTO((int fType));
extern int IDL_CDECL IDL_SetWait IDL_ARG_PROTO((int fType, int iVal));
extern char *IDL_CDECL IDL_MakeTempStruct IDL_ARG_PROTO((IDL_StructDefPtr
        sdef, int n_dim, IDL_MEMINT *dim, IDL_VPTR *var, int zero));
extern char *IDL_CDECL IDL_MakeTempStructVector
        IDL_ARG_PROTO((IDL_StructDefPtr sdef, IDL_MEMINT dim, IDL_VPTR
        *var, int zero));
extern IDL_StructDefPtr IDL_CDECL IDL_MakeStruct IDL_ARG_PROTO((char
        *name, IDL_STRUCT_TAG_DEF *tags));
extern IDL_MEMINT IDL_CDECL IDL_StructTagInfoByName
        IDL_ARG_PROTO((IDL_StructDefPtr sdef, char *name, int msg_action,
        IDL_VPTR *var));
extern IDL_MEMINT IDL_CDECL IDL_StructTagInfoByIndex
        IDL_ARG_PROTO((IDL_StructDefPtr sdef, int index, int msg_action,
        IDL_VPTR *var));
extern char *IDL_CDECL IDL_StructTagNameByIndex
        IDL_ARG_PROTO((IDL_StructDefPtr sdef, int index, int msg_action,
        char **struct_name));
extern int IDL_CDECL IDL_StructNumTags IDL_ARG_PROTO((IDL_StructDefPtr
        sdef));
extern void IDL_CDECL IDL_unform_io IDL_ARG_PROTO((int type, int argc,
        IDL_VPTR *argv, char *argk));
extern void IDL_CDECL IDL_Print IDL_ARG_PROTO((int argc, IDL_VPTR *argv,
        char *argk));
extern void IDL_CDECL IDL_PrintF IDL_ARG_PROTO((int argc, IDL_VPTR *argv,
        char *argk));
extern void IDL_CDECL IDL_Win32MessageLoop  IDL_ARG_PROTO((int fFlush));
extern void IDL_CDECL IDL_RgbToHsv IDL_ARG_PROTO((UCHAR *r, UCHAR *g,
        UCHAR *b, float *h, float *s, float *v, int n));
extern void IDL_CDECL IDL_RgbToHls IDL_ARG_PROTO((UCHAR *r, UCHAR *g,
        UCHAR *b, float *h, float *l, float *s, int n));
extern int IDL_CDECL IDL_AddDevice IDL_ARG_PROTO(( IDL_DEVICE_DEF *dev, 
        int msg_action));
extern char *IDL_OutputFormat[];
extern char *IDL_CDECL IDL_OutputFormatFunc IDL_ARG_PROTO((int type));
extern int IDL_OutputFormatLen[];
extern int IDL_CDECL IDL_OutputFormatLenFunc IDL_ARG_PROTO((int type));
extern char *IDL_OutputFormatNatural[];
extern IDL_LONG IDL_TypeSize[];
extern int IDL_CDECL IDL_TypeSizeFunc IDL_ARG_PROTO((int type));
extern char *IDL_TypeName[];
extern char *IDL_CDECL IDL_TypeNameFunc IDL_ARG_PROTO((int type));
extern IDL_ALLTYPES IDL_zero;
extern IDL_VPTR IDL_CDECL IDL_nonavailable_rtn IDL_ARG_PROTO((int argc,
        IDL_VPTR argv[], char *argk));
extern IDL_SYS_VERSION IDL_SysvVersion;
extern IDL_STRING *IDL_CDECL IDL_SysvVersionArch IDL_ARG_PROTO((void));
extern IDL_STRING *IDL_CDECL IDL_SysvVersionOS IDL_ARG_PROTO((void));
extern IDL_STRING *IDL_CDECL IDL_SysvVersionOSFamily IDL_ARG_PROTO((void));
extern IDL_STRING *IDL_CDECL IDL_SysvVersionRelease IDL_ARG_PROTO((void));
extern char *IDL_ProgramName;
extern char *IDL_CDECL IDL_ProgramNameFunc IDL_ARG_PROTO((void));
extern char *IDL_ProgramNameLC;
extern char *IDL_CDECL IDL_ProgramNameLCFunc IDL_ARG_PROTO((void));
extern IDL_STRING IDL_SysvDir;
extern IDL_STRING *IDL_CDECL IDL_SysvDirFunc IDL_ARG_PROTO((void));
extern IDL_LONG IDL_SysvErrCode;
extern IDL_LONG IDL_CDECL IDL_SysvErrCodeValue IDL_ARG_PROTO((void));
extern IDL_SYS_ERROR_STATE IDL_SysvErrorState;
extern IDL_SYS_ERROR_STATE *IDL_CDECL IDL_SysvErrorStateAddr
        IDL_ARG_PROTO((void));
extern IDL_STRING *IDL_CDECL IDL_SysvErrStringFunc IDL_ARG_PROTO((void));
extern IDL_STRING *IDL_CDECL IDL_SysvSyserrStringFunc
        IDL_ARG_PROTO((void));
extern IDL_LONG IDL_CDECL IDL_SysvErrorCodeValue IDL_ARG_PROTO((void));
extern IDL_LONG *IDL_CDECL IDL_SysvSyserrorCodesAddr IDL_ARG_PROTO((void));
extern IDL_LONG IDL_SysvOrder;
extern IDL_LONG IDL_CDECL IDL_SysvOrderValue IDL_ARG_PROTO((void));
extern int IDL_CDECL IDL_AddSystemRoutine IDL_ARG_PROTO((IDL_SYSFUN_DEF
        *defs, int is_function, int cnt));
extern int IDL_CDECL IDL_BailOut IDL_ARG_PROTO((int stop));
extern int IDL_CDECL IDL_Cleanup IDL_ARG_PROTO((int just_cleanup));
extern int IDL_CDECL IDL_Init IDL_ARG_PROTO((int options, int *argc, char
        *argv[]));
extern int IDL_CDECL IDL_Main IDL_ARG_PROTO((int init_options, int argc,
        char *argv[]));
extern int IDL_CDECL IDL_ExecuteStr IDL_ARG_PROTO((char *cmd));
extern int IDL_CDECL IDL_Execute IDL_ARG_PROTO((int argc, char *argv[]));
extern int IDL_CDECL IDL_RuntimeExec IDL_ARG_PROTO((char *file));
extern void IDL_CDECL IDL_Runtime IDL_ARG_PROTO((int options, int *argc,
        char *argv[], char *file));
extern IDL_TERMINFO IDL_FileTerm;
extern char *IDL_CDECL IDL_FileTermName IDL_ARG_PROTO((void));
extern int IDL_CDECL IDL_FileTermIsTty IDL_ARG_PROTO((void));
extern int IDL_CDECL IDL_FileTermLines IDL_ARG_PROTO((void));
extern int IDL_CDECL IDL_FileTermColumns IDL_ARG_PROTO((void));
extern int IDL_CDECL IDL_FileEnsureStatus IDL_ARG_PROTO((int action, int
        unit, int flags));
extern void IDL_CDECL IDL_FileSetMode IDL_ARG_PROTO((int unit, int
        binary));
extern int IDL_CDECL IDL_FileOpen IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[], char *argk, int access_mode, int extra_flags, int
        longjmp_safe,  int msg_action));
extern void IDL_CDECL IDL_FileClose IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[], char *argk));
extern void IDL_CDECL IDL_FileFlushUnit IDL_ARG_PROTO((int unit));
extern void IDL_CDECL IDL_FileGetUnit IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern void IDL_CDECL IDL_FileFreeUnit IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern int IDL_CDECL IDL_FileSetPtr IDL_ARG_PROTO((int unit, IDL_FILEINT
        pos, int extend, int msg_action));
extern int IDL_CDECL IDL_FileEOF IDL_ARG_PROTO((int unit));
extern void IDL_CDECL IDL_FileStat IDL_ARG_PROTO((int unit, IDL_FILE_STAT
        *stat_blk));
extern void IDL_CDECL IDL_FileSetClose IDL_ARG_PROTO((int unit, int
        allow));
extern IDL_VPTR IDL_CDECL IDL_FileVaxFloat IDL_ARG_PROTO((int argc,
        IDL_VPTR *argv, char *argk));
extern void IDL_CDECL IDL_VarCopy IDL_ARG_PROTO((IDL_REGISTER IDL_VPTR
        src, IDL_REGISTER IDL_VPTR dst));
extern void IDL_CDECL IDL_StoreScalar IDL_ARG_PROTO((IDL_VPTR dest, int
        type, IDL_ALLTYPES *value));
extern void IDL_CDECL IDL_StoreScalarZero IDL_ARG_PROTO((IDL_VPTR dest,
        int type));
extern int IDL_CDECL IDL_KWGetParams IDL_ARG_PROTO((int argc, IDL_VPTR
        *argv, char *argk, IDL_KW_PAR *kw_list, IDL_VPTR *plain_args,  int
        imask));
extern void IDL_CDECL IDL_KWCleanup IDL_ARG_PROTO((int fcn));
extern char *IDL_DitherMethodNames[];
extern char *IDL_CDECL IDL_DitherMethodNamesFunc IDL_ARG_PROTO((int
        method));
extern void IDL_CDECL IDL_RasterDrawThick IDL_ARG_PROTO((IDL_GR_PT *p0,
        IDL_GR_PT *p1, IDL_ATTR_STRUCT *a, IDL_PRO_PTR routine, int
        dot_width));
extern void IDL_CDECL IDL_RasterPolyfill IDL_ARG_PROTO((int *x, int *y,
        int n, IDL_POLYFILL_ATTR *p, IDL_RASTER_DEF *r));
extern void IDL_CDECL IDL_RasterDraw IDL_ARG_PROTO((IDL_GR_PT *p0,
        IDL_GR_PT *p1, IDL_ATTR_STRUCT *a, IDL_RASTER_DEF *r));
extern void IDL_CDECL IDL_Dither IDL_ARG_PROTO((UCHAR *data, int ncols,
        int nrows, IDL_RASTER_DEF *r, int x0, int y0, IDL_TV_STRUCT
        *secondary));
extern void IDL_CDECL IDL_BitmapLandscape IDL_ARG_PROTO((IDL_RASTER_DEF
        *in, IDL_RASTER_DEF *out, int y0));
extern IDL_LONG IDL_CDECL IDL_LongScalar IDL_ARG_PROTO((IDL_REGISTER
        IDL_VPTR p));
extern IDL_ULONG IDL_CDECL IDL_ULongScalar IDL_ARG_PROTO((IDL_REGISTER
        IDL_VPTR p));
extern IDL_LONG64 IDL_CDECL IDL_Long64Scalar IDL_ARG_PROTO((IDL_REGISTER
        IDL_VPTR p));
extern IDL_ULONG64 IDL_CDECL IDL_ULong64Scalar IDL_ARG_PROTO((IDL_REGISTER
        IDL_VPTR p));
extern double IDL_CDECL IDL_DoubleScalar IDL_ARG_PROTO((IDL_REGISTER
        IDL_VPTR p));
extern IDL_MEMINT IDL_CDECL IDL_MEMINTScalar IDL_ARG_PROTO((IDL_REGISTER
        IDL_VPTR p));
extern IDL_FILEINT IDL_CDECL IDL_FILEINTScalar IDL_ARG_PROTO((IDL_REGISTER
        IDL_VPTR p));
extern IDL_VPTR IDL_CDECL IDL_BasicTypeConversion IDL_ARG_PROTO((int argc,
        IDL_VPTR argv[], IDL_REGISTER int type));
extern IDL_VPTR IDL_CDECL IDL_CvtByte IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtFix IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtLng IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtFlt IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtDbl IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtUInt IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtULng IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtLng64 IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtULng64 IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtMEMINT IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtFILEINT IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtComplex IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtDComplex IDL_ARG_PROTO((int argc,
        IDL_VPTR argv[]));
extern IDL_VPTR IDL_CDECL IDL_CvtString IDL_ARG_PROTO((int argc, IDL_VPTR
        argv[], char *argk));
extern int IDL_CDECL IDL_GetKbrd IDL_ARG_PROTO((int should_wait));
extern void IDL_CDECL IDL_VarGetData IDL_ARG_PROTO((IDL_VPTR v, IDL_MEMINT
        *n, char **pd, int ensure_simple));
extern IDL_VPTR IDL_CDECL IDL_ImportArray IDL_ARG_PROTO((int n_dim,
        IDL_MEMINT dim[], int type, UCHAR *data, IDL_ARRAY_FREE_CB
        free_cb,  IDL_StructDefPtr s));
extern IDL_VPTR IDL_CDECL IDL_ImportNamedArray IDL_ARG_PROTO((char *name,
        int n_dim, IDL_MEMINT dim[],  int type, UCHAR *data, 
        IDL_ARRAY_FREE_CB free_cb,  IDL_StructDefPtr s));
extern void IDL_CDECL IDL_Delvar IDL_ARG_PROTO((IDL_VPTR var));
extern void IDL_CDECL IDL_VarEnsureSimple IDL_ARG_PROTO((IDL_VPTR v));
extern IDL_VPTR IDL_CDECL IDL_CvtBytscl IDL_ARG_PROTO((int argc, IDL_VPTR
        *argv, char *argk));
extern void IDL_CDECL IDL_Freetmp IDL_ARG_PROTO((IDL_REGISTER IDL_VPTR p));
extern void IDL_CDECL IDL_Deltmp IDL_ARG_PROTO((IDL_REGISTER IDL_VPTR p));
extern IDL_VPTR IDL_CDECL IDL_Gettmp IDL_ARG_PROTO((void));
extern IDL_VPTR IDL_CDECL IDL_GettmpInt IDL_ARG_PROTO((short value));
extern IDL_VPTR IDL_CDECL IDL_GettmpUInt IDL_ARG_PROTO((IDL_UINT value));
extern IDL_VPTR IDL_CDECL IDL_GettmpLong IDL_ARG_PROTO((IDL_LONG value));
extern IDL_VPTR IDL_CDECL IDL_GettmpULong IDL_ARG_PROTO((IDL_ULONG value));
extern IDL_VPTR IDL_CDECL IDL_GettmpFILEINT IDL_ARG_PROTO((IDL_FILEINT
        value));
extern IDL_VPTR IDL_CDECL IDL_GettmpMEMINT IDL_ARG_PROTO((IDL_MEMINT
        value));
extern char *IDL_CDECL IDL_GetScratch IDL_ARG_PROTO((IDL_REGISTER IDL_VPTR
        *p, IDL_REGISTER IDL_MEMINT n_elts,  IDL_REGISTER IDL_MEMINT
        elt_size));


#ifdef __cplusplus
    }
#endif

#endif                               /* export_IDL_DEF */
