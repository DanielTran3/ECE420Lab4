all: seq.c send_and_recv.c sendrecv.c bcast.c main.c
	mpicc -g -Wall -o seq seq.c Lab4_IO.c -lm
	mpicc -g -Wall -o send_and_recv send_and_recv.c Lab4_IO.c -lm
	mpicc -g -Wall -o sendrecv sendrecv.c Lab4_IO.c -lm
	mpicc -g -Wall -o bcast bcast.c Lab4_IO.c -lm
	mpicc -g -Wall -o main main.c Lab4_IO.c -lm

seq: seq.c
	mpicc -g -Wall -o seq seq.c Lab4_IO.c -lm
send_and_recv: send_and_recv.c
	mpicc -g -Wall -o send_and_recv send_and_recv.c Lab4_IO.c -lm
sendrecv: sendrecv.c
	mpicc -g -Wall -o sendrecv sendrecv.c Lab4_IO.c -lm
bcast: bcast.c
	mpicc -g -Wall -o bcast bcast.c Lab4_IO.c -lm
main: main.c
	mpicc -g -Wall -o main main.c Lab4_IO.c -lm
