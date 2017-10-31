KPL/TM

SPP/WISPR Development Meta-Kernel
===========================================================================

   This meta-kernel loads the necessary SPICE kernels for developing and 
   testing WISPR IDL routines.

Version and Date
---------------------------------------------------------------

   The TEXT_KERNEL_ID stores version information of loaded project text
   kernels.  Each entry associated with the keyword is a string that
   consists of four parts: the kernel name, version, entry date, and type.
   For example, the frame kernel might have an entry as follows:

      \begindata
      TEXT_KERNEL_ID += 'WISPR V0.0.0 16-JUN-2016 TM'
      \begintext

   Version 0.0.0 -- June 16, 2016 -- Angelos Vourlidas

References
---------------------------------------------------------------

      1.   "Kernel Required Reading"

Contact Information
---------------------------------------------------------------

   Direct questions, comments, or concerns about the contents of this kernel
   to:

      Angelos Vourlidas, JHUAPL, (240)228-5073, Angelos.Vourlidas@jhuapl.edu


Implementation Notes
---------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make use
   of this frame kernel must "load" the kernel normally during program
   initialization.  Loading the kernel associates the data items with
   their names in a data structure called the "kernel pool".  The SPICELIB
   routine FURNSH loads a kernel into the pool as shown below:

      FORTRAN: (SPICELIB)

         CALL FURNSH ( frame_kernel_name )

      C: (CSPICE)

         furnsh_c ( frame_kernel_name );

      IDL: (ICY)

         cspice_furnsh, frame_kernel_name

   In order for a program or routine to extract data from the pool, the
   SPICELIB routines GDPOOL, GIPOOL, and GCPOOL are used.  See [2] for
   more details.

   This file was created and may be updated with a text editor or word
   processor.


Kernels
---------------------------------------------------------------

   The following kernels are loaded with this kernel file:

      Kernel Name		    	Type     
      ============= 		   ================      

      de421.bsp			   Solar Ephemeris                   
      pck00010.tpc      	   Planetary Constants
      naif0011.tls		   Leapseconds
      heliospheric.tf		   Dynamic Heliospheric Frames for the NASA STEREO mission
      SPP20180731P2_ephem.bsp	   Latest SPP s/c kernel
      spp_v004.tf		   SPP dynamic frame 
      spp_wispr_v001.ti		   WISPR Instrument Kernel

   \begindata
	  PATH_VALUES = ('/usr/local/ssw/packages/sunspice/data',
		         '/data1/idl/kernels/gen',
			 '/data1/idl/solohi/kernels/SPP',
			 '/data1/idl/solohi/kernels/SOLO')

      PATH_SYMBOLS = ( 'GEN', 'SITE', 'SPP', 'SOLO')

      KERNELS_TO_LOAD = ( '$GEN/pck00010.tpc',
                          '$GEN/naif0011.tls',
                          '$GEN/de405.bsp',  
      			  '$GEN/heliospheric.tf',
      		          '$SPP/SPP20180731P2_ephem.bsp',
                          '$SPP/spp_v004.tf',
			  '$SPP/spp_rtn.tf',
      		          '$SPP/spp_wispr_v002.ti'               ) 
   \begintext



      

    