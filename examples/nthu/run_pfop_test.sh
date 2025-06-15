#!/bin/bash

# All paths relative to executable
./pchtrees --ntrees 20 --mphalo 1e12  --params ./pfop_test.toml --zmax 30.0 --no-output-trees --process-first-order-progenitors
