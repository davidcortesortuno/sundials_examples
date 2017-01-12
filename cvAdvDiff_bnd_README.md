This example was compiled by doing:

1. Compile the libraries:

   ```
   gcc -c src/sundials/sundials_math.c -Iinclude 
   gcc -c src/sundials/sundials_band.c -Iinclude 
   ```

   This generates the `.o` files (the -c flag compiles but does not link)

2. Compile the example file:

   ```
   gcc -Wall -c cvAdvDiff_bnd.c -o test.o -Iinclude -L$LD_LIBRARY_PATH -lm -lsundials_cvodes -lsundials_ida -lsundials_kinsol -lsundials_nvecserial
   ```

   This generates the `test.o` file

3. Now link all the `.o` files including the libraries and headers:

   ```
   gcc -Wall test.o sundials_math.o sundials_band.o -o test -Iinclude -L$LD_LIBRARY_PATH -lm -lsundials_cvodes -lsundials_ida -lsundials_kinsol -lsundials_nvecserial
   ```

4. Execute the `test`: `./test` 


## Notes


We can create the shared libraries:

```
gcc -shared -o libsundials_math.so sundials_math.o
gcc -shared -o libsundials_band.so sundials_band.o
```

And then use the libraries to directly compile the example C file:

```
gcc -Wall cvAdvDiff_bnd.c -o test -L/home/david/tmp/cvode -lsundials_math -lsundials_band -L$LD_LIBRARY_PATH -lm -lsundials_cvodes -lsundials_ida -lsundials_kinsol -lsundials_nvecserial -Iinclude
```

Somehow, order matters. The `sundials_*` libraries must be called before the
ones from `LD_LIBRARY_PATH`
