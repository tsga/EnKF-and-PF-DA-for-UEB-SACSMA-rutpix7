MAINTAINERCLEANFILES	= Makefile.in

AM_CXXFLAGS =  -Wall	-I$(top_builddir)	\
			-I$(BOOSTINCLUDEDIR)	\
			-I$(top_builddir)/HydroDomain	\
			-I$(top_builddir)/ModelBase		\
			-I$(top_builddir)/Utilities		\
			-I$(top_builddir)/Interface		\
			-I$(top_builddir)/Input		\
			-I$(top_builddir)/Output	\
			-I$(top_builddir)/hlrms		\
			-I$(top_srcdir)		\
			-I$(top_builddir)/loki                \
                        -I../HydroDomain                \
                        -I../ModelBase                  \
                        -I../Utilities                  \
                        -I../Interface                  \
                        -I../Input                      \
                        -I../Output                     \
                        -I../hlrms                      \
                        -I.

SOURCES 	=	sacFuncInside_Ens.cpp		\
			uebDAsrPFdafunctions.cpp	\
                        uebDAsrPFFuncBeforeLoop.cpp	\
                       	uebDAsrPFFuncInsideLoop.cpp     

CXX_SRCS =  	fland1_Ens.cpp			\
		sacCudaCalls.cpp		\
		uebDAsrPFcudacalls.cpp		\
		uebDAsrPFuebdecls.cpp			\
		uebDAsrPFsnowdv.cpp			\
		uebDAsrPFsnowdgtv.cpp			\
		uebDAsrPFcanopy.cpp			\
		uebDAsrPFsnowxv.cpp			\
		uebDAsrPFuebinputs.cpp			\
		uebDAsrPFncfunctions.cpp		\
		uebDAsrPFmatrixnvector.cpp

OBJS = $(CXX_SRCS:.cpp=.o)

DLINOBJ = uebDAPFGPU.o
		
						
UEBDALIBCPPFLAGS 	= -I$(NETCDFC_INCDIR) -I$(CUDA_INCDIR) -I$(BOOSTINCLUDEDIR) -I$HOME/codes/include -std=c++11 
UEBDALIBLDFLAGS 	= -L$(NETCDFC_LIBDIR) -L$(CUDA_LIBDIR) -L$(HOME)/codes/lib -L. $(DLINOBJ) $(OBJS) -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -ldl -lcuda -lcudart -lcudadevrt
                        
lib_LIBRARIES 		= libuebDAsrPF.a
libuebDAsrPF_a_SOURCES 	= $(SOURCES)
libuebDAsrPF_a_CPPFLAGS = $(UEBDALIBCPPFLAGS)
libuebDAsrPF_a_LIBADD = $(DLINOBJ) $(OBJS)