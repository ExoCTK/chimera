# Environment Setup

To run the retrieval code on the multi-core science servers, `ssh` into a machine (with names science1 through science7) with

```
ssh -XY <username>@science6
```

Next, make sure you are in your root development environment with

```
source deactivate
```

Then `cd` into the Linux directory containing the `pyMultiNest` install with

```
cd /grp/software/Linux/multinest_test_python2.7
```

Finally, add the `MultiNest` library path and the modules path to your environment with

```
export LD_LIBRARY_PATH="/grp/software/Linux/multinest_test/anaconda3/lib"
export PATH=/grp/software/Linux/multinest_test_python2.7/anaconda2/bin:$PATH
export PYTHONPATH=<path_to_local_chimera_package>:$PYTHONPATH
```
And you're all set!

## Compile!
```
rm _tran_module.so
python setup.py build_ext —inplace
rm cea2.x
./compile_cea_64bit.com
```
