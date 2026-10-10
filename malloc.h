#ifndef MALLOC_H
#define MALLOC_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#ifdef __linux__
#include <sys/mman.h>
#endif

#define HUGEPAGE ((size_t)1 << 21)

static void *_hugealloc(size_t alignment, size_t size) {
  void *ret;
  if(size > SIZE_MAX - (HUGEPAGE - 1))
    return NULL;
  size = (size + HUGEPAGE - 1) & ~(HUGEPAGE - 1);
  if(posix_memalign(&ret,(alignment > HUGEPAGE) ? alignment : HUGEPAGE,size))
    return NULL;
#ifdef MADV_HUGEPAGE
  madvise(ret,size,MADV_HUGEPAGE);
#endif
  return ret;
}

static void *_malloc(size_t size) {
  void *ret;
  ret = (size >= HUGEPAGE) ? _hugealloc(0,size) : malloc(size);
  if(ret == NULL) {
    fprintf(stderr,"ERROR: Not enough memory\n");
    exit(1);
  }
  return ret;
}

static void *_aligned_alloc(size_t alignment, size_t size) {
  void *ret;
  ret = (size >= HUGEPAGE) ? _hugealloc(alignment,size) : aligned_alloc(alignment,size);
  if(ret == NULL) {
    fprintf(stderr,"ERROR: Not enough memory\n");
    exit(1);
  }
  return ret;
}

static void *_aligned_calloc(size_t alignment, size_t size) {
  void *ret;
  ret = (size >= HUGEPAGE) ? _hugealloc(alignment,size) : aligned_alloc(alignment,size);
  if(ret == NULL) {
    fprintf(stderr,"ERROR: Not enough memory\n");
    exit(1);
  }
  memset(ret,0,size);
  return ret;
}

static void *_realloc(void *ptr, size_t size) {
  void *ret;
  ret = realloc(ptr, size);
  if(ret == NULL) {
    fprintf(stderr,"ERROR: Not enough memory\n");
    exit(1);
  }
  return ret;
}

#endif
