#Compiler config for the main target
CC = g++ -std=c++11
LD = $(CC) -std=c++11
CFLAGS = -Wall -Wextra -O2 -I /usr/local/include
#CFLAGS += -Wall -Wextra -Werror -pedantic -ansi -Wshadow -Wdouble-promotion -Wundef -fno-common -Wconversion -Wunused-parameter
TEST_CFLAGS += $(CFLAGS) -I$(FUSED_GTEST_TMP_DIR) -larmadillo -Og -DGTEST_HAS_PTHREAD=0
LDFLAGS = -Wall -Wextra -larmadillo

#Modules to consider in the build. foo.cpp will be foo.
include tests/modules
include src/modules

#Folders config
BINDIR = bin
OBJDIR = obj
SRCDIR = src
DOCDIR = doc
TEST_SRCDIR = tests
FUSED_GTEST_TMP_DIR = tmp
GTEST_SRC = gtest

#Names of the targets
TARGET = $(BINDIR)/solver
TEST_TARGET = $(BINDIR)/tests

all : makedirs $(TARGET)

makedirs :
	mkdir -p tmp
	mkdir -p bin
	mkdir -p obj

#Build of the target
MAIN_SRC = $(addprefix $(SRCDIR)/, $(MAIN:=.cpp))
MAIN_OBJ = $(addprefix $(OBJDIR)/, $(MAIN:=.o))
SOURCES = $(addprefix $(SRCDIR)/, $(MODULES:=.cpp))
OBJECTS = $(addprefix $(OBJDIR)/, $(MODULES:=.o))
ALL_OBJECTS = $(OBJECTS) $(MAIN_OBJ)
ALL_HEADERS = $(addprefix $(SRCDIR)/, $(ORPHANED_HEADERS:=.h))
ALL_SOURCES = $(SOURCES) $(MAIN_SRC)

$(TARGET) : $(ALL_OBJECTS) $(ALL_SOURCES) $(ALL_HEADERS)
	$(LD) $(LDFLAGS) -o $@ $(ALL_OBJECTS)

$(MAIN_OBJ): obj/%.o : src/%.cpp $(MAIN_SRC)
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJECTS): obj/%.o : src/%.cpp $(SOURCES) $(ALL_HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<


#Doxygen
.PHONY : doc
doc :
	doxygen $(DOCDIR)/Doxyfile


FUSED_GTEST_H = $(FUSED_GTEST_TMP_DIR)/gtest/gtest.h
FUSED_GTEST_ALL_CC = $(FUSED_GTEST_TMP_DIR)/gtest/gtest-all.cc
GTEST_MAIN_CC = $(GTEST_SRC)/googletest/src/gtest_main.cc
TEST_SOURCES = $(addprefix $(TEST_SRCDIR)/, $(TEST_MODULES:=.cpp))
TEST_OBJECTS = $(addprefix $(OBJDIR)/, $(TEST_MODULES:=.o))
ALL_TEST_OBJECTS = $(OBJECTS) $(OBJDIR)/gtest-all.o $(OBJDIR)/gtest_main.o $(TEST_OBJECTS)

.PHONY : tests
tests : makedirs $(TEST_TARGET)

check : makedirs $(TEST_TARGET)
	$(TEST_TARGET) --gtest_color=yes

#GTEST special objects
$(FUSED_GTEST_H) :
	$(GTEST_SRC)/googletest/scripts/fuse_gtest_files.py $(FUSED_GTEST_TMP_DIR)

$(FUSED_GTEST_ALL_CC) :
	$(GTEST_SRC)/googletest/scripts/fuse_gtest_files.py $(FUSED_GTEST_TMP_DIR)

$(OBJDIR)/gtest-all.o : $(FUSED_GTEST_H) $(FUSED_GTEST_ALL_CC)
	$(CC) $(TEST_CFLAGS) -c $(FUSED_GTEST_TMP_DIR)/gtest/gtest-all.cc -o $(OBJDIR)/gtest-all.o

$(OBJDIR)/gtest_main.o : $(FUSED_GTEST_H) $(GTEST_MAIN_CC)
	$(CC) $(TEST_CFLAGS) -c $(GTEST_MAIN_CC) -o $(OBJDIR)/gtest_main.o

$(TEST_OBJECTS): $(OBJDIR)%.o : $(TEST_SRCDIR)%.cpp $(TEST_SOURCES)
	$(CC) $(TEST_CFLAGS) -c -o $@ $< -I $(GTEST_SRC)/googletest/include

$(TEST_TARGET) : $(ALL_TEST_OBJECTS)
	$(LD) $(TEST_CFLAGS) $^ -o $(TEST_TARGET)

.PHONY : clean
clean :
	rm -rf $(ALL_TEST_OBJECTS) $(ALL_OBJECTS)
	rm -rf $(TARGET) $(TEST_TARGET)
	rm -rf $(FUSED_GTEST_TMP_DIR)
	rm -rf $(DOCDIR)/html

coefs:
	python3 scripts/generateHermiteCoefs.py 100 > src/hermiteCoefs.h
