htmldir = @htmldir@/html
htmlsearchdir = $(htmldir)/search

HTML_INST = html/*.css html/*.html html/*.png 
SEARCH_INST = html/search/*

if ENABLE_DOXYGEN
  html_DATA = $(HTML_INST)
  htmlsearch_DATA = $(SEARCH_INST)
else
  html_DATA =
  htmlsearch_DATA = 
endif

.PHONY: $(HTML_INST) $(SEARCH_INST)

$(HTML_INST): html/index.html
doc-local: html-local

html/index.html: meshkit.dox
	$(DOXYGEN) $<

uninstall-local:
	cd $(DESTDIR) && rmdir $(htmlsearchdir) && rmdir $(htmldir) || :

