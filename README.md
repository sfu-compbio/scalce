## SCALCE - Sequence Compression Algorithms using Locally Consistent Encoding

### How to install?

In most cases, typing

	git clone https://github.com/sfu-compbio/scalce.git
	cd scalce
	make download
	make

should result in a happy `scalce` executable.
`make download` will download all necessary dependencies necessary for building SCALCE.

We also recommend that you install [pigz](http://zlib.net/pigz/) and make the `pigz`
executable reachable via `PATH` variable.

If this doesn't work, please open a Github issue and we'll try to sort it out.

### Warning

Variable-length reads (e.g. Pacbio) are not fully supported! There is a limited support (described
in documentation); however, it is not tested and is rather buggy. Use caution and 
make sure to test decompression validity on Pacbio reads!


### More stuff...

... is available at http://sfu-compbio.github.io/scalce/.
