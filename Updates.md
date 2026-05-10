Updates:
gfortran --version
gfortran -O3 -ffixed-line-length-none ram.f -o ram.exe

Once ram is built, can run ./ram.exe < ram.in or whatever input file. Right now, ram.in is hardcoded as the input file it is looking for, so you need to:cp ram1.in ram.in
./ram.exe . Where ram1.in is whatever your chosen input file is. 
