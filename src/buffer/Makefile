all:
	gcc -c -Wall -Werror -fPIC buffer.c -o buffer.o
	gcc -c -Wall -Werror -fPIC buffer2d.c -o buffer2d.o
	gcc -shared buffer.o buffer2d.o -o libbuffer.so
	gcc -c -Wall -Werror test.c -o test.o
	gcc test.o -L. -lbuffer -o test
