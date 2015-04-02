VERSION := 2.0.0
DIR = asmvar

align.so:

	cd $(DIR); gcc align.c -fPIC -shared -o align.so

clean:
	
	rm $(DIR)/align.so

.PHONY: clean
