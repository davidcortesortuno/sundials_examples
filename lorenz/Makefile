#!/bin/bash
# Assuming sundials is installed:
all:
	make compile
	make run
	make plot
compile:
	gcc -Wall lorenz.c -o lorenz -lsundials_cvodes -lsundials_nvecserial -lm
run:
	./lorenz > lorenz.txt
plot:
	python plot_system.py
