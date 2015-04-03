#include "Stack.h"

void push(item** stack_p, double d, int index)
{
    item* tmp = (*stack_p);
    (*stack_p) = (item*)malloc(sizeof(item));
    (*stack_p)->val = d;
    (*stack_p)->ndx = index;
    (*stack_p)->next = tmp;
}

int pop(item** stack_p, thread_data* p)
{
    item* tmp = (*stack_p)->next;
    if ((*stack_p) == NULL) return -1;
    p->val = (*stack_p)->val;
    p->ndx = (*stack_p)->ndx;
    free(*stack_p);
    (*stack_p) = tmp;
    return 0;
}

