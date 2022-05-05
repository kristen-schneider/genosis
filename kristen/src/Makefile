CXX=g++
LIBS=-lhts
RM=rm -f

SRCS=main.cpp \
     utils.cpp \
     slice.cpp \
     read_encodings.cpp \
     faiss_pm.cpp

OBJS=$(subst .cpp,.o,$(SRCS))

all: k

k: $(OBJS)
	$(CXX) -o k $(OBJS) $(LIBS)

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) k
