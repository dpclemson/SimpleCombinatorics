RFLAGS = `root-config --cflags --libs`

CFLAGS = -Wall

PROGRAMS = do_rangenx

all: $(PROGRAMS)

clean:
	rm -f $(PROGRAMS)

debug: CFLAGS += -g
debug: $(PROGRAMS)

phenix: CFLAGS += -std=c++11
phenix: $(PROGRAMS)

do_rangenx: do_rangenx.C
	g++ -o do_rangenx do_rangenx.C $(CFLAGS) $(RFLAGS)

