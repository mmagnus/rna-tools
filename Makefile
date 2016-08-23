.PHONY: gh-pages
gh-pages:
	git checkout gh-pages
	rm -rf build _sources _static _modules source
	git checkout master $(GH_PAGES_SOURCES)
	git reset HEAD
	make -f docs/Makefile html
	rm -rf $(GH_PAGES_SOURCES) build
	git add -A
	git commit -m "Generated gh-pages for `git log master -1 --pretty=short \
    --abbrev-commit`" && git push origin gh-pages ; git checkout master
