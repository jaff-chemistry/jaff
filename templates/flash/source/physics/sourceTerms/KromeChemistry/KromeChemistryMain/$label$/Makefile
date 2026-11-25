# ABOUTME: Makefile of the JAFF-generated network for the KromeChemistry unit; FLASH's
# ABOUTME: ODEPACK (opkd*.o) comes from the parent KromeChemistryMain Makefile

KromeChemistry += user_functions.o chemistry.o chemistry_dlsodes.o \
                  KromeChemistry_data.o KromeChemistry_init.o KromeChemistry.o

chemistry.o           : user_functions.o
chemistry_dlsodes.o   : chemistry.o
KromeChemistry_init.o : KromeChemistry_data.o chemistry.o
KromeChemistry.o      : KromeChemistry_data.o chemistry.o chemistry_dlsodes.o \
	Timers_interface.o Grid_interface.o Eos_interface.o
