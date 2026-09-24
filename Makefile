CXX = g++
SOURCES = src/FileUtils.cpp barcodeNinja.cpp
TARGET = barcodeNinja
LIBDIR = lib

BACKEND ?= C
OBJDIR = build/$(BACKEND)
OBJECTS = $(SOURCES:%.cpp=$(OBJDIR)/%.o)

ifeq ($(BACKEND),R)

    R_HOME := $(shell R RHOME)
    R_VERSION := $(shell R --version | head -n 1 | awk '{print $$3}')

    CPPFLAGS := -std=c++17 -O2 -Wall -Wextra -DRCPP_ACTIVE=0 -DLIBR_ACTIVE=1 \
                -DR_HOME_COMPILED=\"$(R_HOME)\" \
                -DR_VERSION_COMPILED=\"$(R_VERSION)\" \
                $(shell R CMD config --cppflags) \
                -I$(R_HOME)/include

    LDFLAGS := $(shell R CMD config --ldflags) \
               -Wl,-rpath,'$$ORIGIN/$(LIBDIR)' \
               -lz

else ifeq ($(BACKEND),C)

    CPPFLAGS := -std=c++17 -O2 -Wall -Wextra -DRCPP_ACTIVE=0 -DLIBR_ACTIVE=0
    LDFLAGS  := -lz

else

    $(error BACKEND must be C or R, got '$(BACKEND)')

endif

$(TARGET): $(OBJECTS)
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)
ifeq ($(BACKEND),R)
	$(MAKE) bundle-libs
	$(MAKE) bundle-rhome
endif

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) -MMD -MP -c $< -o $@

-include $(OBJECTS:.o=.d)

bundle-libs:
	@echo "Bundling libR and dependencies..."
	@mkdir -p $(LIBDIR)
	cp -v $(R_HOME)/lib/libR.so $(LIBDIR)
	ldd $(TARGET) | awk '{print $$3}' | grep -E '^/' | \
	grep -vE '^/lib|^/usr/lib|ld-linux' | \
	while read lib; do \
	    echo "Copying $$lib"; \
	    cp -n $$lib $(LIBDIR) 2>/dev/null || true; \
	done

bundle-rhome:
	@echo "Bundling full R runtime..."
	@mkdir -p $(LIBDIR)
	rsync -a $(R_HOME)/ $(LIBDIR)/R/
	@test -f $(LIBDIR)/R/etc/Renviron || (echo "ERROR: Missing Renviron" && exit 1)
	@test -d $(LIBDIR)/R/library/base || (echo "ERROR: Missing base package" && exit 1)
	@echo "R runtime successfully bundled"

clean:
	rm -rf build $(LIBDIR) $(TARGET)
	rm -f barcodeNinja.o src/FileUtils.o

.PHONY: clean bundle-libs bundle-rhome
