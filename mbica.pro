
TEMPLATE = subdirs
SUBDIRS = library \
		app \
                tests

library.subdir = source
app.depends = library
tests.depends = library

