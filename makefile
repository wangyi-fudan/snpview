snpview: snpview.cpp
	g++ snpview.cpp -o snpview -O3 -Wall -lbam -Lsamtools-0.1.19/ -Isamtools-0.1.19/ -static -s -lpthread -lz

