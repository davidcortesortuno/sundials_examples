#!/bin/bash
# Assuming sundials is installed:
all:
	make compile
	make run
	make plot
compile:
	gcc -Wall lorenz.c -o lorenz -O3 -I/usr/include -L/usr/lib -lsundials_cvode -lsundials_nvecserial -lm
run:
	./lorenz > lorenz.txt
plot:
	python plot_system.py
