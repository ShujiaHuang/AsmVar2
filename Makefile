VERSION := 2.0.0
DIR = asmvar

all: align.so encode.so

align.so:

	cd $(DIR); gcc align.c -fPIC -shared -o align.so

encode.so:

	cd $(DIR)/utils; g++ encode.cpp -fPIC -shared -o encode.so

clean:
	
	rm $(DIR)/align.so $(DIR)/utils/encode.so

.PHONY: clean
