
MFEM_DIR = ../../contrib/mfem-4.6

include $(MFEM_DIR)/config/config.mk
EIGEN_INCLUDE := -I../../contrib/eigen-3.4.0

#gnu compiler
CXX		= $(MFEM_CXX)

#Ubuntu
CXXFLAGS	= $(MFEM_FLAGS)

INCLUDES	= $(wildcard ../../src/*.h)
SRCS		= $(wildcard ../../src/*.cpp)
OBJS		= $(patsubst %.cpp, %.o, $(SRCS))

DC_INCLUDE	= -I../../src/

%.o : %.cpp
	@echo "IGAFEMDC is compiling "$<"..."
	@$(CXX) $(CXXFLAGS) $(DC_INCLUDE) $(MFEM_INCLUDE) $(EIGEN_INCLUDE) -c $< -o $@ 

IGAFEMDC: $(OBJS)
	@$(CXX) $(CXXFLAGS) -o IGAFEMDC $(OBJS) $(MFEM_LIBS) 


clean:
	@rm -rf $(OBJS)
	
cleanall: 
	@rm -rf $(OBJS) IGAFEMDC
