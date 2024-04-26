#include <stdlib.h>
#include <string.h>
#include <assert.h>
 
#include "list.h"

list_t *list_new(size_t elementSize, freeFunction_t freeFunc, const char *name)
{
  assert(elementSize > 0);

  list_t *list = (list_t*) malloc(sizeof(list_t));

  list->logicalLength = 0;
  list->elementSize = elementSize;
  list->head = list->tail = NULL;
  list->freeFunc = freeFunc;
  list->name = name ? strdup(name) : NULL;

  return list;
}

void list_destroy(list_t *list)
{
  if ( list )
    {
      listNode_t *current;
      while ( list->head )
        {
          current = list->head;
          list->head = current->next;
        
          if ( list->freeFunc )
            list->freeFunc(current->data);

          free(current->data);
          free(current);
        }
      if ( list->name ) free(list->name);
      free(list);
    }
}

void list_prepend(list_t *list, void *element)
{
  listNode_t *node = (listNode_t *) malloc(sizeof(listNode_t));
  node->data = malloc(list->elementSize);
  memcpy(node->data, element, list->elementSize);

  node->next = list->head;
  list->head = node;

  if ( !list->tail ) list->tail = list->head;

  list->logicalLength++;
}

void list_append(list_t *list, void *element)
{
  listNode_t *node = (listNode_t *) malloc(sizeof(listNode_t));
  node->data = malloc(list->elementSize);
  node->next = NULL;

  memcpy(node->data, element, list->elementSize);

  if ( list->logicalLength == 0 )
    {
      list->head = list->tail = node;
    }
  else
    {
      list->tail->next = node;
      list->tail = node;
    }

  list->logicalLength++;
}

void list_for_each(list_t *list, listIterator_t iterator)
{
  assert(iterator != NULL);
 
  listNode_t *node = list->head;
  bool result = true;
  while ( node && result )
    {
      result = iterator(node->data);
      node = node->next;
    }
}

void list_head(list_t *list, void *element, bool removeFromList)
{
  assert(list->head != NULL);
 
  listNode_t *node = list->head;
  memcpy(element, node->data, list->elementSize);
 
  if ( removeFromList )
    {
      list->head = node->next;
      list->logicalLength--;
 
      free(node->data);
      free(node);
    }
}

void list_tail(list_t *list, void *element)
{
  assert(list->tail != NULL);
  listNode_t *node = list->tail;
  memcpy(element, node->data, list->elementSize);
}

void *list_entry(list_t *list, int index)
{
  if ( list )
    {
      int i = 0;
      listNode_t *node = list->head;
      while ( node )
        {
          if ( i == index ) return *(void **)node->data;
          node = node->next;
          ++i;
        }
    }

  return NULL;
}

int list_size(list_t *list)
{
  return list->logicalLength;
}

const char *list_name(list_t *list)
{
  return list->name;
}
