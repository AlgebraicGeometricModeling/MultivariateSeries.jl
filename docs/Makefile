.SUFFIXES: .ipynb .md

MD := $(patsubst %.ipynb,%.md,$(wildcard mrkd/expl/*.ipynb))

.ipynb.md:
	jupyter nbconvert --to markdown $*.ipynb
	perl -pi -e "s/\\$$\\$$/\\$$/g" $*.md

all: $(MD) MultivariateSeries.jl
	julia --color=yes make.jl

MultivariateSeries.jl:
	git clone -b gh-pages https://github.com/AlgebraicGeometricModeling/MultivariateSeries.jl.git

publish:
	cd MultivariateSeries.jl; git add docs; git commit -am "update gh-pages"; git push

clean:
	rm ./mrkd/expl/*.md



