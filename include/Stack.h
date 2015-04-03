#ifndef STACK_H
#define STACK_H
#include <stdlib.h>
#include <string.h>


typedef struct elem 
{
    double val;
    int ndx;
    struct elem* next;
} item;

typedef struct data
{
    double val;
    int ndx;
} thread_data;


void push(item** stack_p, double d, int index);
int pop(item** stack_p, thread_data* p);

#endif
