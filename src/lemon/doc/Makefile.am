EXTRA_DIST = \
	Doxyfile.in \
	DoxygenLayout.xml \
	coding_style.dox \
	dirs.dox \
	groups.dox \
	lgf.dox \
	license.dox \
	mainpage.dox \
	migration.dox \
	min_cost_flow.dox \
	named-param.dox \
	namespaces.dox \
	html \
	CMakeLists.txt

DOC_EPS_IMAGES18 = \
	grid_graph.eps \
	nodeshape_0.eps \
	nodeshape_1.eps \
	nodeshape_2.eps \
	nodeshape_3.eps \
	nodeshape_4.eps

DOC_EPS_IMAGES27 = \
	bipartite_matching.eps \
	bipartite_partitions.eps \
	connected_components.eps \
	edge_biconnected_components.eps \
	matching.eps \
	node_biconnected_components.eps \
	planar.eps \
	strongly_connected_components.eps

DOC_EPS_IMAGES = \
	$(DOC_EPS_IMAGES18) \
	$(DOC_EPS_IMAGES27)

DOC_PNG_IMAGES = \
	$(DOC_EPS_IMAGES:%.eps=gen-images/%.png)

EXTRA_DIST += $(DOC_EPS_IMAGES:%=images/%)

#html:
#	$(MAKE) $(AM_MAKEFLAGS) html

GS_COMMAND=gs -dNOPAUSE -dBATCH -q -dEPSCrop -dTextAlphaBits=4 -dGraphicsAlphaBits=4

$(DOC_EPS_IMAGES18:%.eps=gen-images/%.png): gen-images/%.png: images/%.eps
	-mkdir gen-images
	$(GS_COMMAND) -sDEVICE=pngalpha -r18 -sOutputFile=$@ $<;

$(DOC_EPS_IMAGES27:%.eps=gen-images/%.png): gen-images/%.png: images/%.eps
	-mkdir gen-images
	$(GS_COMMAND) -sDEVICE=pngalpha -r27 -sOutputFile=$@ $<;

references.dox: references.bib
	  python @abs_top_srcdir@/scripts/bib2dox.py @abs_top_builddir@/$< >$@;

html-local: $(DOC_PNG_IMAGES) 
	  doxygen Doxyfile;

clean-local:
	-rm -rf html
	-rm -f doxygen.log
	-rm -f $(DOC_PNG_IMAGES)
	-rm -rf gen-images

update-external-tags:
	wget -O libstdc++.tag.tmp http://gcc.gnu.org/onlinedocs/libstdc++/latest-doxygen/libstdc++.tag && \
	mv libstdc++.tag.tmp libstdc++.tag || \
	rm libstdc++.tag.tmp

install-html-local: html
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(htmldir)/html
	for p in html/*.{html,css,png,map,gif,tag} ; do \
	  f="`echo $$p | sed -e 's|^.*/||'`"; \
	  echo " $(INSTALL_DATA) $$p $(DESTDIR)$(htmldir)/html/$$f"; \
	  $(INSTALL_DATA) $$p $(DESTDIR)$(htmldir)/html/$$f; \
	done

uninstall-local:
	@$(NORMAL_UNINSTALL)
	for p in html/*.{html,css,png,map,gif,tag} ; do \
	  f="`echo $$p | sed -e 's|^.*/||'`"; \
	  echo " rm -f $(DESTDIR)$(htmldir)/html/$$f"; \
	  rm -f $(DESTDIR)$(htmldir)/html/$$f; \
	done

.PHONY: update-external-tags
