The directory "Complex" contains a Makefile and code for testing complex data type with OpenACC. It defines the gauge fields as complex doubles. PGI does not like this and it will not compile. The compiler throws the same errors as with the 2D-Schwinger model code that I attempted to accelerate

The directory "Real" is the same code, but defines the gauge fields as doubles. This compiles and runs. 
