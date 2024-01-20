cd /home/musab/Desktop/tsp_openmp

gcc tsp_openmp.c -fopenmp

{ time ./a.out 1 < 11.in >> ./output/tsp_1t_11.txt;  } 2>> ./output/tsp_1t_11.txt;
{ time ./a.out 1 < 14.in >> ./output/tsp_1t_14.txt;  } 2>> ./output/tsp_1t_14.txt;
{ time ./a.out 1 < 17.in >> ./output/tsp_1t_17.txt;  } 2>> ./output/tsp_1t_17.txt;
{ time ./a.out 1 < 20.in >> ./output/tsp_1t_20.txt;  } 2>> ./output/tsp_1t_20.txt;
{ time ./a.out 1 < 23.in >> ./output/tsp_1t_23.txt;  } 2>> ./output/tsp_1t_23.txt;

{ time ./a.out 4 < 11.in >> ./output/tsp_4t_11.txt;  } 2>> ./output/tsp_4t_11.txt;
{ time ./a.out 4 < 14.in >> ./output/tsp_4t_14.txt;  } 2>> ./output/tsp_4t_14.txt;
{ time ./a.out 4 < 17.in >> ./output/tsp_4t_17.txt;  } 2>> ./output/tsp_4t_17.txt;
{ time ./a.out 4 < 20.in >> ./output/tsp_4t_20.txt;  } 2>> ./output/tsp_4t_20.txt;
{ time ./a.out 4 < 23.in >> ./output/tsp_4t_23.txt;  } 2>> ./output/tsp_4t_23.txt;

{ time ./a.out 8 < 11.in >> ./output/tsp_8t_11.txt;  } 2>> ./output/tsp_8t_11.txt;
{ time ./a.out 8 < 14.in >> ./output/tsp_8t_14.txt;  } 2>> ./output/tsp_8t_14.txt;
{ time ./a.out 8 < 17.in >> ./output/tsp_8t_17.txt;  } 2>> ./output/tsp_8t_17.txt;
{ time ./a.out 8 < 20.in >> ./output/tsp_8t_20.txt;  } 2>> ./output/tsp_8t_20.txt;
{ time ./a.out 8 < 23.in >> ./output/tsp_8t_23.txt;  } 2>> ./output/tsp_8t_23.txt;

{ time ./a.out 16 < 11.in >> ./output/tsp_16t_11.txt;  } 2>> ./output/tsp_16t_11.txt;
{ time ./a.out 16 < 14.in >> ./output/tsp_16t_14.txt;  } 2>> ./output/tsp_16t_14.txt;
{ time ./a.out 16 < 17.in >> ./output/tsp_16t_17.txt;  } 2>> ./output/tsp_16t_17.txt;
{ time ./a.out 16 < 20.in >> ./output/tsp_16t_20.txt;  } 2>> ./output/tsp_16t_20.txt;
{ time ./a.out 16 < 23.in >> ./output/tsp_16t_23.txt;  } 2>> ./output/tsp_16t_23.txt;
