---
title: "Notes"
---

# Release notes

- *(10-Jan-2016)* **SCALCE version 2.8 release**
	- Bugfixes (arithmetic decoding bugfix)
	- Fixed a decompression bug when number of reads was greater than 2^32. Compression was not affected.
	- *New*: support for variable length reads via `scalce-pacbio`.


- *(20-May-2013)* **SCALCE version 2.7 release**
	- Bugfixes (no-arithmetic fix)

- *(13-May-2013)* **SCALCE version 2.6 release**
	- Bugfixes

- *(02-Apr-2013)* **SCALCE version 2.5 release**
	- Read splitting supported
	- Standard output during decompression supported
	- Bugfixes

- *(25-Mar-2013)* **SCALCE version 2.4 release**
	- Auto-pigz detection
	- Bugfixes

- *(10-Sep-2012)* **SCALCE version 2.3 release**
	- Decompression speed improvements

- *(25-Jul-2012)* **SCALCE version 2.2 release**
	- Speed improvements
	- Arithmetic coding for qualities is now optional
	- Multiple bug fixes

- *(06-Jun-2012)* **SCALCE version 2.1 release**
	- Better compression of reads
	- Arithmetic coding for qualities
	- Multiple bug fixes

- *(02-Mar-2012)* **SCALCE version 1.4 release**
	- Serious data loss when using multithreading bug fixed

- *(20-Feb-2012)* **SCALCE version 1.3 release**
	- Various bug fixes

- *(17-Feb-2012)* **SCALCE version 1.2 release**
	- Various bug fixes

- *(08-Feb-2012)* **SCALCE version 1.1 release**
	- OpenMP support
	- pigz support

- *(06-Dec-2011)* **SCALCE version 1.0 release**
	- Initial release
