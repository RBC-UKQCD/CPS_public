SUBDIRS=rpcgenvml 

.PHONY: $(SUBDIRS)

all: $(SUBDIRS)

clean: $(SUBDIRS)

$(SUBDIRS):
	echo recursing into $(SUBDIRS)
	$(MAKE) -C $@ $(MAKECMDGOALS)
