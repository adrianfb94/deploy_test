#ifndef _LIST_H
#define _LIST_H

#include <stdio.h>
#include <stdbool.h>

// a common function used to free malloc'd objects
typedef void (*freeFunction_t)(void *);
 
typedef bool (*listIterator_t)(void *);
 
typedef struct _listNode_t {
  void *data;
  struct _listNode_t *next;
} listNode_t;
 
typedef struct {
  size_t logicalLength;
  size_t elementSize;
  char *name;
  listNode_t *head;
  listNode_t *tail;
  freeFunction_t freeFunc;
} list_t;
 
list_t *list_new(size_t elementSize, freeFunction_t freeFunc, const char *name);
void list_destroy(list_t *list);
 
void list_prepend(list_t *list, void *element);
void list_append(list_t *list, void *element);
int list_size(list_t *list);
const char *list_name(list_t *list);

void list_for_each(list_t *list, listIterator_t iterator);
void list_head(list_t *list, void *element, bool removeFromList);
void list_tail(list_t *list, void *element);
void *list_entry(list_t *list, int index);

#endif
