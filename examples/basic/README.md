To run the test in this directory, first compile a pchtrees executable (see the
main README). Then run `./run_basic_example.sh` in your shell.

The ./pchtree symlink in this examples directory will only work if your main
directory is two levels above this one, and contains your pchtree binary. This
should be the case if you're working in a github checkout and you complied
without problems. If the example doesn't work, either remake the symlink to
point to your binary, or copy your binary to this directory.

The same applies to the `./data` directory, which is also a symlink in this
example. Note that the path to the data directory can also be set in the
parameter file (the example here assumes the path is `./data`)
