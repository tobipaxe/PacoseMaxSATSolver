all: preprocessor

preprocessor:
	$(MAKE) -C src all
	mv src/maxpre maxpre

.PHONY: lib
lib:
	$(MAKE) -C src lib

clean:
	$(MAKE) -C src clean
