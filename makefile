all:
	upcc -gupc -network=udp -pthreads=4 rooms.upc -o rooms

source:
	source /opt/nfs/config/source_bupc.sh

run:
	UPC_NODEFILE=nodes upcrun -n 32 ./rooms 10

clean:
	rm rooms nodes

